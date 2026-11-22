// Joint k/epsilon Newton direction for the ORIGINAL physical HB residual.
// Relative field coordinates improve scaling; prescribed epsilon wall rows
// are eliminated and reconstructed from k on EVERY residual evaluation.
// U and face flux remain fixed. No pseudo-time term enters this Jacobian.
{
    labelList offsets(nO+1,0);
    scalar volume=0;
    forAll(subTimeMeshes,K)
    {
        offsets[K+1]=offsets[K]+2*subTimeMeshes[K].nCells();
        for (const scalar v : subTimeMeshes[K].V()) volume+=v;
    }
    reduce(volume,sumOp<scalar>());
    const label size=offsets[nO];
    scalarField weights(size,Zero),inverseDiagonal(size,Zero);
    PtrList<scalarField> spatialLower(2*nO),spatialUpper(2*nO);
    const auto finite = [&](const scalarField& x)
    {
        bool ok=true;
        for (const scalar value : x) ok=ok && std::isfinite(value);
        reduce(ok,andOp<bool>()); return ok;
    };
    const auto dotProduct = [&](const scalarField& a,const scalarField& b)
    {
        scalar value=0;
        forAll(a,i) value+=weights[i]*a[i]*b[i];
        reduce(value,sumOp<scalar>()); return value;
    };
    const auto norm = [&](const scalarField& x)
    {
        const scalar v=dotProduct(x,x);
        return std::isfinite(v) && v>=0 ? Foam::sqrt(v) : scalar(GREAT);
    };
    forAll(subTimeMeshes,K)
    {
        const HBPhaseTime clock(runTime,snapshots[K]+hbMotionDeltaT,
            hbMotionDeltaT,subTimeMeshes[K],meshPhiHB[K]);
        for (label f=0;f<2;++f)
        {
            const bool eps=f==1;
            tmp<fvScalarMatrix> equation=hbTurbulence[K].phaseImplicitEquation(eps);
            hbTurbulence[K].phaseImplicitConstrain(equation.ref(),eps);
            spatialLower.set(2*K+f,new scalarField(subTimeMeshes[K].nInternalFaces(),Zero));
            spatialUpper.set(2*K+f,new scalarField(subTimeMeshes[K].nInternalFaces(),Zero));
            if (equation().hasLower() || equation().hasUpper())
            {
                spatialLower[2*K+f]=equation().lower();
                spatialUpper[2*K+f]=equation().upper();
            }
            const volScalarField diagonal(equation().A());
            const scalar scale=eps ? phaseEpsilonScale : phaseKScale;
            const scalarField& base=eps ? phaseTransactionBase[K].epsilon.primitiveField()
                : phaseTransactionBase[K].k.primitiveField();
            forAll(base,i)
            {
                const label index=offsets[K]+2*i+f;
                weights[index]=subTimeMeshes[K].V()[i]/volume;
                if (eps && phaseTransactionEpsilonWall[K][i]) continue;
                if (!std::isfinite(diagonal[i]) || diagonal[i]<=0)
                    FatalErrorInFunction << "Nonpositive turbulence Newton preconditioner diagonal"
                        << exit(FatalError);
                inverseDiagonal[index]=scale/(diagonal[i]*base[i]);
            }
        }
    }
    phaseRestoreTurbulenceBase();
    // Retain the stiff HB phase coupling in the right preconditioner.
    // These are LOCAL cell inverses, not global spatial matrix inverses.
    const label cells=subTimeMeshes[0].nCells();
    List<scalarField> phaseInverse(2*cells);
    {
        PtrList<scalarField> vd(nO*nO);
        forAll(subTimeMeshes,K)
            forAll(subTimeMeshes,J)
            {
                blockFvMatrix<scalar,scalar> block(subTimeMeshes[K]);
                block.diag()=Zero; HB.addBlock(block,K,J);
                vd.set(K*nO+J,new scalarField(block.diag()));
            }
        scalar worstError=0;
        for (label cell=0;cell<cells;++cell)
            for (label f=0;f<2;++f)
            {
                scalarSquareMatrix matrix(nO,Zero),inverse(nO,Zero);
                const bool eps=f==1;
                const scalar scale=eps ? phaseEpsilonScale : phaseKScale;
                forAll(subTimeMeshes,K)
                {
                    if (eps && phaseTransactionEpsilonWall[K][cell])
                    { matrix(K,K)=1; continue; }
                    matrix(K,K)=1/inverseDiagonal[offsets[K]+2*cell+f];
                    forAll(subTimeMeshes,J)
                    {
                        if (eps && phaseTransactionEpsilonWall[J][cell]) continue;
                        const scalar base=eps ? phaseTransactionBase[J].epsilon[cell]
                            : phaseTransactionBase[J].k[cell];
                        matrix(K,J)+=vd[K*nO+J][cell]*base
                            /(subTimeMeshes[K].V()[cell]*scale);
                    }
                }
                LUscalarMatrix(matrix).inv(inverse);
                scalarField& saved=phaseInverse[2*cell+f];
                saved.setSize(nO*nO);
                scalar aNorm=0,iNorm=0,errorNorm=0;
                for (label i=0;i<nO;++i)
                {
                    scalar aSum=0,iSum=0,eSum=0;
                    for (label j=0;j<nO;++j)
                    {
                        saved[i*nO+j]=inverse(i,j);
                        aSum+=mag(matrix(i,j)); iSum+=mag(inverse(i,j));
                        scalar product=0;
                        for (label k=0;k<nO;++k) product+=matrix(i,k)*inverse(k,j);
                        eSum+=mag(product-(i==j ? 1.0 : 0.0));
                    }
                    aNorm=max(aNorm,aSum); iNorm=max(iNorm,iSum);
                    errorNorm=max(errorNorm,eSum);
                }
                const scalar error=errorNorm/max(aNorm*iNorm,scalar(VSMALL));
                if (!std::isfinite(error) || error>1e-10)
                    FatalErrorInFunction << "Invalid local turbulence HB inverse"
                        << exit(FatalError);
                worstError=max(worstError,error);
            }
        reduce(worstError,maxOp<scalar>());
        Info<< "IHBJointTurbulencePreconditioner phaseBlockInverseBackwardError "
            << worstError << " phaseBlockSize " << nO << endl;
    }
    // Symmetric block Gauss-Seidel on native INTERNAL spatial faces. The
    // diagonal blocks retain HB; processor/AMI actions remain in the true J.
    // Topology correspondence was checked before entering this update.
    const labelUList& owner=phaseReferenceMesh.lduAddr().lowerAddr();
    const labelUList& neighbour=phaseReferenceMesh.lduAddr().upperAddr();
    labelList forwardStart(cells+1,0),backwardStart(cells+1,0);
    forAll(owner,face)
    {
        ++forwardStart[min(owner[face],neighbour[face])+1];
        ++backwardStart[max(owner[face],neighbour[face])+1];
    }
    for (label cell=0;cell<cells;++cell)
    { forwardStart[cell+1]+=forwardStart[cell]; backwardStart[cell+1]+=backwardStart[cell]; }
    labelList forwardFaces(owner.size()),backwardFaces(owner.size());
    labelList forwardNext(forwardStart),backwardNext(backwardStart);
    forAll(owner,face)
    {
        forwardFaces[forwardNext[min(owner[face],neighbour[face])]++]=face;
        backwardFaces[backwardNext[max(owner[face],neighbour[face])]++]=face;
    }
    const auto precondition = [&](const scalarField& input,scalarField& output)
    {
        output=scalarField(size,Zero);
        scalarField work(input);
        const auto solveCell = [&](const label cell)
        {
            for (label f=0;f<2;++f)
            {
                const scalarField& inverse=phaseInverse[2*cell+f];
                forAll(subTimeMeshes,K)
                {
                    if (f==1 && phaseTransactionEpsilonWall[K][cell]) continue;
                    scalar value=0;
                    forAll(subTimeMeshes,J)
                        value+=inverse[K*nO+J]*work[offsets[J]+2*cell+f];
                    output[offsets[K]+2*cell+f]=value;
                }
            }
        };
        const auto distribute = [&](const label cell,const bool forward)
        {
            const labelList& start=forward ? forwardStart : backwardStart;
            const labelList& faces=forward ? forwardFaces : backwardFaces;
            for (label position=start[cell];position<start[cell+1];++position)
            {
                const label face=faces[position];
                const bool currentOwner=owner[face]==cell;
                const label other=currentOwner ? neighbour[face] : owner[face];
                forAll(subTimeMeshes,K)
                    for (label f=0;f<2;++f)
                    {
                        const scalar coefficient=currentOwner ? spatialLower[2*K+f][face]
                            : spatialUpper[2*K+f][face];
                        const scalar base=f==1 ? phaseTransactionBase[K].epsilon[cell]
                            : phaseTransactionBase[K].k[cell];
                        const scalar scale=f==1 ? phaseEpsilonScale : phaseKScale;
                        work[offsets[K]+2*other+f]-=coefficient*base
                            *output[offsets[K]+2*cell+f]/(subTimeMeshes[K].V()[other]*scale);
                    }
            }
        };
        // After each forward solve work_i=D_i*y_i, providing the middle D
        // in (D+U)^-1 D (D+L)^-1 without an extra matrix multiplication.
        for (label cell=0;cell<cells;++cell)
        { solveCell(cell); distribute(cell,true); }
        for (label cell=cells-1;cell>=0;--cell)
        { solveCell(cell); distribute(cell,false); }
    };
    const auto packResidual = [&]()
    {
        scalarField result(size,Zero);
        forAll(subTimeMeshes,K)
        {
            for (label f=0;f<2;++f)
            {
                const bool eps=f==1;
                const scalar scale=eps ? phaseEpsilonScale : phaseKScale;
                const scalarField& residual=eps ? phaseNewtonEpsilonResidual[K]
                    : phaseNewtonKResidual[K];
                forAll(residual,i)
                    if (!(eps && phaseTransactionEpsilonWall[K][i]))
                        result[offsets[K]+2*i+f]=residual[i]/scale;
            }
        }
        return result;
    };
    scalar kr=0,er=0;
    phaseProposalValid=phaseMeasureTurbulence(kr,er);
    const scalarField rhs(packResidual());
    const scalar initial=norm(rhs);
    const scalar target=phaseTurbulenceRelTol*initial;
    label evaluations=0;
    const auto applyJacobian = [&](const scalarField& x,scalarField& result)
    {
        scalar maximum=0;
        for (const scalar value : x) maximum=max(maximum,mag(value));
        reduce(maximum,maxOp<scalar>());
        if (maximum<=VSMALL) { result=scalarField(size,Zero); return true; }
        const scalar step=1e-6/maximum;
        phaseRestoreTurbulenceBase();
        forAll(subTimeMeshes,K)
            for (label f=0;f<2;++f)
            {
                const bool eps=f==1;
                scalarField& q=hbTurbulence[K].phaseImplicitField(eps).primitiveFieldRef();
                const scalarField& base=eps ? phaseTransactionBase[K].epsilon.primitiveField()
                    : phaseTransactionBase[K].k.primitiveField();
                forAll(q,i)
                    if (!(eps && phaseTransactionEpsilonWall[K][i]))
                        q[i]=base[i]*(1+step*x[offsets[K]+2*i+f]);
            }
        scalar k=0,e=0;
        const bool valid=phaseMeasureTurbulence(k,e);
        ++evaluations;
        if (valid) result=(rhs-packResidual())/step;
        phaseRestoreTurbulenceBase();
        return valid && finite(result);
    };

    scalarField increment(size,Zero),residual(rhs);
    scalar trueRelative=1;
    label directionsUsed=0;
    Info<< "IHBJointTurbulenceNewtonStart initialScaledRms " << initial
        << " targetScaledRms " << target << " fieldCoordinates relative"
        << " physicalJacobian 1 eliminatedEpsilonWallRows 1" << endl;
    for (label cycle=0;cycle<phaseTurbulenceCycles && phaseProposalValid;++cycle)
    {
        const scalar beta=norm(residual);
        if (beta<=target) break;
        const label m=phaseTurbulenceDirections;
        List<scalarField> basis(m+1),preconditioned(m);
        List<scalarField> h(m+1,scalarField(m,Zero));
        scalarField cosine(m,Zero),sine(m,Zero),g(m+1,Zero);
        basis[0]=residual/beta; g[0]=beta;
        label used=0;
        for (label j=0;j<m;++j)
        {
            precondition(basis[j],preconditioned[j]);
            scalarField action(size,Zero);
            if (!applyJacobian(preconditioned[j],action))
            { phaseProposalValid=false; break; }
            // Reorthogonalize because the FD operator has finite accuracy.
            for (label pass=0;pass<2;++pass)
                for (label i=0;i<=j;++i)
                {
                    const scalar a=dotProduct(basis[i],action);
                    h[i][j]+=a; action-=a*basis[i];
                }
            h[j+1][j]=norm(action);
            const bool breakdown=h[j+1][j]<=VSMALL;
            if (!breakdown) basis[j+1]=action/h[j+1][j];
            for (label i=0;i<j;++i)
            {
                const scalar a=cosine[i]*h[i][j]+sine[i]*h[i+1][j];
                h[i+1][j]=-sine[i]*h[i][j]+cosine[i]*h[i+1][j];
                h[i][j]=a;
            }
            const scalar divisor=std::hypot(h[j][j],h[j+1][j]);
            if (!std::isfinite(divisor) || divisor<=VSMALL)
            { phaseProposalValid=false; break; }
            cosine[j]=h[j][j]/divisor; sine[j]=h[j+1][j]/divisor;
            h[j][j]=divisor; h[j+1][j]=0;
            g[j+1]=-sine[j]*g[j]; g[j]*=cosine[j]; used=j+1;
            if (mag(g[j+1])<=target || breakdown) break;
        }
        if (!phaseProposalValid || !used) break;
        scalarField y(used,Zero);
        for (label i=used-1;i>=0;--i)
        {
            scalar value=g[i];
            for (label j=i+1;j<used;++j) value-=h[i][j]*y[j];
            y[i]=value/h[i][i];
        }
        for (label i=0;i<used;++i) increment+=y[i]*preconditioned[i];
        scalarField action(size,Zero);
        if (!finite(increment) || !applyJacobian(increment,action))
        { phaseProposalValid=false; break; }
        residual=rhs-action; trueRelative=norm(residual)/max(initial,scalar(VSMALL));
        directionsUsed+=used;
        Info<< "IHBJointTurbulenceKrylov cycle " << cycle+1
            << " directions " << used << " trueRelative " << trueRelative
            << " targetRelative " << phaseTurbulenceRelTol << endl;
    }
    phaseRestoreTurbulenceBase();
    scalar damping=1;
    if (phaseProposalValid)
    {
        forAll(subTimeMeshes,K)
            for (label f=0;f<2;++f)
            {
                const bool eps=f==1;
                const scalar minimum=hbTurbulence[K].phaseImplicitMinimum(eps);
                const scalarField& base=eps ? phaseTransactionBase[K].epsilon.primitiveField()
                    : phaseTransactionBase[K].k.primitiveField();
                forAll(base,i)
                {
                    if (eps && phaseTransactionEpsilonWall[K][i]) continue;
                    const scalar dq=base[i]*increment[offsets[K]+2*i+f];
                    if (dq<0) damping=min(damping,0.99*(base[i]-minimum)/(-dq));
                }
            }
        reduce(damping,minOp<scalar>());
        phaseProposalValid=std::isfinite(damping) && damping>=phaseTurbulenceMinDamping;
    }
    if (phaseProposalValid)
        forAll(subTimeMeshes,K)
            for (label f=0;f<2;++f)
            {
                const bool eps=f==1;
                scalarField& q=hbTurbulence[K].phaseImplicitField(eps).primitiveFieldRef();
                const scalarField& base=eps ? phaseTransactionBase[K].epsilon.primitiveField()
                    : phaseTransactionBase[K].k.primitiveField();
                forAll(q,i)
                    if (!(eps && phaseTransactionEpsilonWall[K][i]))
                        q[i]=base[i]*(1+damping*increment[offsets[K]+2*i+f]);
            }
    Info<< "IHBJointTurbulenceNewtonEnd finite " << phaseProposalValid
        << " trueRelative " << trueRelative << " linearConverged "
        << (trueRelative<=phaseTurbulenceRelTol) << " directions " << directionsUsed
        << " evaluations " << evaluations << " positivityDamping " << damping
        << " proposalOnly 1" << endl;
}
