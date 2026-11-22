// Application-local pressure Schur RIGHT preconditioner for stage FGMRES.
// Included after nkScalarWork/nkVectorWork (physical integrated RHS) exist.
// M=[S C; G A]. Use a full cell-local velocity/HB inverse B, and approximate
// S-C*B*G by native -laplacian(rAUphysical,pCorrection), including the gauge.
// The actual J, H, and face-flux discretisation are untouched.
List<scalarField> nkSchurVelocityInverse;
PtrList<volScalarField> nkSchurPressure(nO);
PtrList<volVectorField> nkSchurVelocity(nO);
PtrList<fvScalarMatrix> nkSchurPressureMatrix(nO);
PtrList<fvScalarMatrix::fvSolver> nkSchurPressureSolver(nO);
PtrList<vectorField> nkSchurGradient(nO);
PtrList<labelList> nkSchurForwardStart(nO),nkSchurForwardFaces(nO);
PtrList<labelList> nkSchurBackwardStart(nO),nkSchurBackwardFaces(nO);
label nkSchurSweepChecks=0;
dictionary nkSchurControls;
const label nkSchurVelocitySize=3*nO;
bool nkSchurIsotropic=true;

if (stagePreconditioner=="pressureSchur")
{
    const label cellCount=subTimeMeshes[0].nCells();
    forAll(subTimeMeshes,K)
        if (subTimeMeshes[K].nCells()!=cellCount)
            FatalErrorInFunction << "Local HB Schur inverse requires matching phase cells"
                << exit(FatalError);
    // Absolute Cartesian HB and scalar momentum transport usually give
    // A_KJ=a_KJ*I. Its three identical phase blocks need only ONE inverse.
    // Detect exact isotropy; retain the general tensor path for other cases.
    forAll(subTimeMeshes,K)
        forAll(subTimeMeshes,J)
            if (eqSystem.dVByVExists(K,J))
                for (const tensor& a : eqSystem.dVByV(K,J).diag())
                    for (direction cmpt=0;cmpt<9;++cmpt)
                        nkSchurIsotropic=nkSchurIsotropic
                            && a.component(cmpt)==(cmpt%4==0 ? a.xx() : scalar(0));
    reduce(nkSchurIsotropic,andOp<bool>());
    nkSchurVelocityInverse.setSize(cellCount);
    const label n=nkSchurIsotropic ? nO : nkSchurVelocitySize;
    scalar inverseBackwardError=0;
    for (label celli=0;celli<cellCount;++celli)
    {
        scalarField lu(n*n,Zero),rowScale(n,Zero);
        labelList permutation(n);
        forAll(subTimeMeshes,K)
            forAll(subTimeMeshes,J)
                if (eqSystem.dVByVExists(K,J))
                {
                    const tensor& block=eqSystem.dVByV(K,J).diag()[celli];
                    if (nkSchurIsotropic) lu[K*n+J]=block.xx();
                    else
                    for (direction i=0;i<3;++i)
                        for (direction j=0;j<3;++j)
                            lu[(3*K+i)*n+3*J+j]=block.component(3*i+j);
                }
        const scalarField originalVelocityBlock(lu);
        // Invert A itself, never the velocity sub-block of a coupled inverse.
        // Row scaling and partial pivoting follow the validated local LU path.
        for (label row=0;row<n;++row)
        {
            permutation[row]=row;
            scalar largest=0;
            for (label col=0;col<n;++col)
            {
                const scalar value=lu[row*n+col];
                if (!std::isfinite(value))
                    FatalErrorInFunction << "Nonfinite velocity block in cell " << celli
                        << exit(FatalError);
                largest=max(largest,mag(value));
            }
            if (largest<=VSMALL)
                FatalErrorInFunction << "Zero velocity block row in cell " << celli
                    << exit(FatalError);
            rowScale[row]=1.0/largest;
            for (label col=0;col<n;++col) lu[row*n+col]*=rowScale[row];
        }
        for (label col=0;col<n;++col)
        {
            label pivot=col;
            scalar largest=mag(lu[col*n+col]);
            for (label row=col+1;row<n;++row)
                if (mag(lu[row*n+col])>largest)
                {
                    largest=mag(lu[row*n+col]);
                    pivot=row;
                }
            if (!std::isfinite(largest) || largest<=stageSchurPivotTolerance)
                FatalErrorInFunction << "Near-singular Schur velocity block in cell "
                    << celli << ", column " << col << ", scaled pivot " << largest
                    << exit(FatalError);
            if (pivot!=col)
            {
                for (label j=0;j<n;++j)
                {
                    const scalar value=lu[col*n+j];
                    lu[col*n+j]=lu[pivot*n+j];
                    lu[pivot*n+j]=value;
                }
                const label index=permutation[col];
                permutation[col]=permutation[pivot];
                permutation[pivot]=index;
            }
            for (label row=col+1;row<n;++row)
            {
                lu[row*n+col]/=lu[col*n+col];
                for (label j=col+1;j<n;++j)
                    lu[row*n+j]-=lu[row*n+col]*lu[col*n+j];
            }
        }
        scalarField& inverse=nkSchurVelocityInverse[celli];
        inverse.setSize(n*n);
        scalarField solution(n,Zero);
        for (label column=0;column<n;++column)
        {
            for (label row=0;row<n;++row)
            {
                scalar value=permutation[row]==column ? rowScale[column] : 0;
                for (label j=0;j<row;++j) value-=lu[row*n+j]*solution[j];
                solution[row]=value;
            }
            for (label row=n-1;row>=0;--row)
            {
                scalar value=solution[row];
                for (label j=row+1;j<n;++j) value-=lu[row*n+j]*solution[j];
                solution[row]=value/lu[row*n+row];
                if (!std::isfinite(solution[row]))
                    FatalErrorInFunction << "Nonfinite Schur velocity inverse in cell "
                        << celli << exit(FatalError);
                inverse[row*n+column]=solution[row];
            }
        }
        scalar normA=0,normInverse=0,normDefect=0;
        for (label row=0;row<n;++row)
        {
            scalar aSum=0,iSum=0,eSum=0;
            for (label col=0;col<n;++col)
            {
                aSum+=mag(originalVelocityBlock[row*n+col]);
                iSum+=mag(inverse[row*n+col]);
                scalar product=0;
                for (label k=0;k<n;++k)
                    product+=originalVelocityBlock[row*n+k]*inverse[k*n+col];
                eSum+=mag(product-(row==col ? 1.0 : 0.0));
            }
            normA=max(normA,aSum);normInverse=max(normInverse,iSum);
            normDefect=max(normDefect,eSum);
        }
        const scalar relative=normDefect/max(normA*normInverse,scalar(VSMALL));
        if (!std::isfinite(relative) || relative>1e-10)
            FatalErrorInFunction << "Schur velocity inverse check failed in cell "
                << celli << ": " << relative << exit(FatalError);
        inverseBackwardError=max(inverseBackwardError,relative);
    }

    reduce(inverseBackwardError,maxOp<scalar>());
    Info<< "IHBStageSchurInverse backwardError " << inverseBackwardError << endl;

    nkSchurControls.set("solver",word("GAMG"));
    nkSchurControls.set("smoother",word("GaussSeidel"));
    nkSchurControls.set("agglomerator",word("faceAreaPair"));
    nkSchurControls.set("cacheAgglomeration",true);
    nkSchurControls.set("nCellsInCoarsestLevel",label(10));
    nkSchurControls.set("mergeLevels",label(1));
    nkSchurControls.set("nPreSweeps",label(1));
    nkSchurControls.set("nPostSweeps",label(1));
    nkSchurControls.set("nFinestSweeps",label(2));
    // v2112 direct coarse LU does not support cyclicAMI interfaces.
    // Native iterative coarse solves retain AMI; stored Z permits variability.
    nkSchurControls.set("directSolveCoarsest",false);
    nkSchurControls.set("tolerance",scalar(0));
    nkSchurControls.set("relTol",scalar(0));
    nkSchurControls.set("minIter",stageSchurVCycles);
    nkSchurControls.set("maxIter",stageSchurVCycles);
    nkSchurControls.set("log",label(0));

    forAll(subTimeMeshes,K)
    {
        fvMesh& mesh=subTimeMeshes[K];
        const HBPhaseTime phaseClock
        (
            runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,
            mesh,meshPhiHB[K]
        );
        wordList pressureTypes(p[K].boundaryField().types());
        forAll(pressureTypes,patchi)
        {
            const fvPatchScalarField& patch=p[K].boundaryField()[patchi];
            if (!patch.coupled() && patch.type()!="empty"
             && patch.type()!="fixedValue" && patch.type()!="zeroGradient")
                FatalErrorInFunction
                    << "pressureSchur only supports homogeneous fixedValue, zeroGradient, "
                    << "empty and coupled pressure increments; phase " << K
                    << ", patch " << patch.patch().name() << " has type " << patch.type()
                    << exit(FatalError);
        }
        // Existing pIncr has zeroGradient physical patches; it must NOT be used
        // for this native Poisson solve because the fixed-p outlet is essential.
        nkSchurPressure.set(K,new volScalarField
        (
            IOobject("pStageSchur",runTime.timeName(),mesh,
                IOobject::NO_READ,IOobject::NO_WRITE,false),
            mesh,dimensionedScalar("zero",p[K].dimensions(),0),pressureTypes
        ));
        nkSchurVelocity.set(K,new volVectorField
        (
            IOobject("uStageSchur",runTime.timeName(),mesh,
                IOobject::NO_READ,IOobject::NO_WRITE,false),
            mesh,dimensionedVector("zero",U[K].dimensions(),Zero),
            vectorSolIncrement[K].boundaryField().types()
        ));
        nkSchurGradient.set(K,new vectorField(mesh.nCells(),Zero));
        nkSchurPressure[K].correctBoundaryConditions();
        surfaceScalarField response(fvc::interpolate(flowFrozenResponse[K]));
        if (flowFacePseudoTime) response *= flowFaceWeight[K];
        // Select the physical p discretisation even though the scratch field has
        // a unique name. Cache only its implicit compact Laplacian: explicit
        // nonorthogonal terms are not re-injected as a source during P actions.
        const word laplacianName
            ("laplacian("+response.name()+","+p[K].name()+")");
        nkSchurPressureMatrix.set(K,new fvScalarMatrix
        (
            -fvm::laplacian(response,nkSchurPressure[K],laplacianName)
        ));
        fvScalarMatrix& pressureMatrix=nkSchurPressureMatrix[K];
        pressureMatrix.source()=Zero;
        if (pressureNeedsReference[K] && pRefCells[K]>=0)
            pressureMatrix.diag()[pRefCells[K]]+=pressureGaugeWeight[K];
        autoPtr<fvScalarMatrix::fvSolver> pressureSolver
            (pressureMatrix.solver(nkSchurControls));
        nkSchurPressureSolver.set(K,pressureSolver.ptr());

        if (stageSchurVelocitySweeps)
        {
            const blockFvMatrix<vector,tensor>& a=eqSystem.dVByV(K,K);
            const labelUList& owner=a.lduAddr().lowerAddr();
            const labelUList& neighbour=a.lduAddr().upperAddr();
            const label nFaces=(a.hasLower() || a.hasUpper()) ? owner.size() : 0;
            if (nFaces && (a.lower().size()!=nFaces || a.upper().size()!=nFaces
                || neighbour.size()!=nFaces))
                FatalErrorInFunction << "Schur velocity face/address size mismatch"
                    << exit(FatalError);
            nkSchurForwardStart.set(K,new labelList(cellCount+1,0));
            nkSchurBackwardStart.set(K,new labelList(cellCount+1,0));
            labelList& forwardStart=nkSchurForwardStart[K];
            labelList& backwardStart=nkSchurBackwardStart[K];
            for (label facei=0;facei<nFaces;++facei)
            {
                const label lo=min(owner[facei],neighbour[facei]);
                const label hi=max(owner[facei],neighbour[facei]);
                if (lo<0 || hi>=cellCount || lo==hi)
                    FatalErrorInFunction << "Invalid Schur internal face addressing"
                        << exit(FatalError);
                ++forwardStart[lo+1];
                ++backwardStart[hi+1];
            }
            for (label celli=0;celli<cellCount;++celli)
            {
                forwardStart[celli+1]+=forwardStart[celli];
                backwardStart[celli+1]+=backwardStart[celli];
            }
            nkSchurForwardFaces.set(K,new labelList(nFaces));
            nkSchurBackwardFaces.set(K,new labelList(nFaces));
            labelList forwardNext(forwardStart),backwardNext(backwardStart);
            for (label facei=0;facei<nFaces;++facei)
            {
                const label lo=min(owner[facei],neighbour[facei]);
                const label hi=max(owner[facei],neighbour[facei]);
                nkSchurForwardFaces[K][forwardNext[lo]++]=facei;
                nkSchurBackwardFaces[K][backwardNext[hi]++]=facei;
            }
        }
    }
    Info<< "IHBStagePressureSchur phases " << nO
        << " localVelocityBlock " << nkSchurVelocitySize
        << " storedInverseSize " << (nkSchurIsotropic ? nO : nkSchurVelocitySize)
        << " pressureVCycles " << stageSchurVCycles
        << " velocitySweeps " << stageSchurVelocitySweeps
        << " response physical gauge original side right method FGMRES" << endl;
}

const auto nkSchurCellInverse = [&](const label celli,
    const scalarField& rhs,scalarField& result)
{
    const scalarField& inverse=nkSchurVelocityInverse[celli];
    bool finite=true;
    for (label row=0;row<nkSchurVelocitySize;++row)
    {
        scalar value=0;
        if (nkSchurIsotropic)
            for (label phase=0;phase<nO;++phase)
                value+=inverse[(row/3)*nO+phase]*rhs[3*phase+row%3];
        else
            for (label col=0;col<nkSchurVelocitySize;++col)
                value+=inverse[row*nkSchurVelocitySize+col]*rhs[col];
        result[row]=value;
        finite=finite && std::isfinite(value);
    }
    return finite;
};

const auto nkSchurApplyLocalVelocityInverse = [&](PtrList<vectorField>& values)
{
    const label n=nkSchurVelocitySize;
    scalarField rhs(n,Zero),result(n,Zero);
    bool finite=true;
    forAll(nkSchurVelocityInverse,celli)
    {
        forAll(values,K)
            for (direction cmpt=0;cmpt<3;++cmpt)
                rhs[3*K+cmpt]=values[K][celli].component(cmpt);
        const bool cellFinite=nkSchurCellInverse(celli,rhs,result);
        finite=finite && cellFinite;
        forAll(values,K)
            for (direction cmpt=0;cmpt<3;++cmpt)
                values[K][celli].component(cmpt)=result[3*K+cmpt];
    }
    reduce(finite,andOp<bool>());
    return finite;
};

// B=(D+U)^-1 D (D+L)^-1, with D the full 3*nO HB/pseudo-time block.
// Interface coefficients remain in the outer true J; these sweeps cover only
// internal faces. Both Schur velocity actions use exactly this same B.
const auto nkSchurApplyVelocityInverse = [&](PtrList<vectorField>& values)
{
    if (!stageSchurVelocitySweeps) return nkSchurApplyLocalVelocityInverse(values);
    const label n=nkSchurVelocitySize;
    const label cells=nkSchurVelocityInverse.size();
    const bool check=stageSchurCheckSweeps && nkSchurSweepChecks<2;
    PtrList<vectorField> originalRhs(nO),forwardSolution(nO);
    if (check)
        forAll(values,K)
        {
            originalRhs.set(K,new vectorField(values[K]));
            forwardSolution.set(K,new vectorField(cells,Zero));
        }
    scalarField rhs(n,Zero),result(n,Zero);
    bool finite=true;
    const auto solveCell = [&](const label celli)
    {
        forAll(values,K)
            for (direction cmpt=0;cmpt<3;++cmpt)
                rhs[3*K+cmpt]=values[K][celli].component(cmpt);
        const bool cellFinite=nkSchurCellInverse(celli,rhs,result);
        finite=finite && cellFinite;
    };
    // Keep work_i=D_i*y_i after the forward solve; only distribute L*y_i.
    // This supplies the required middle D without another matrix product.
    for (label celli=0;celli<cells;++celli)
    {
        solveCell(celli);
        forAll(values,K)
        {
            const vector y(result[3*K],result[3*K+1],result[3*K+2]);
            if (check) forwardSolution[K][celli]=y;
            const blockFvMatrix<vector,tensor>& a=eqSystem.dVByV(K,K);
            const labelUList& owner=a.lduAddr().lowerAddr();
            const labelUList& neighbour=a.lduAddr().upperAddr();
            const labelList& start=nkSchurForwardStart[K];
            for (label pos=start[celli];pos<start[celli+1];++pos)
            {
                const label facei=nkSchurForwardFaces[K][pos];
                const bool ownerFirst=owner[facei]<neighbour[facei];
                const label later=ownerFirst ? neighbour[facei] : owner[facei];
                const tensor& coefficient=ownerFirst ? a.lower()[facei] : a.upper()[facei];
                values[K][later]-=dot(coefficient,y);
            }
        }
    }
    for (label celli=cells-1;celli>=0;--celli)
    {
        solveCell(celli);
        forAll(values,K)
        {
            const vector z(result[3*K],result[3*K+1],result[3*K+2]);
            values[K][celli]=z;
            const blockFvMatrix<vector,tensor>& a=eqSystem.dVByV(K,K);
            const labelUList& owner=a.lduAddr().lowerAddr();
            const labelUList& neighbour=a.lduAddr().upperAddr();
            const labelList& start=nkSchurBackwardStart[K];
            for (label pos=start[celli];pos<start[celli+1];++pos)
            {
                const label facei=nkSchurBackwardFaces[K][pos];
                const bool ownerFirst=owner[facei]<neighbour[facei];
                const label earlier=ownerFirst ? owner[facei] : neighbour[facei];
                const tensor& coefficient=ownerFirst ? a.upper()[facei] : a.lower()[facei];
                values[K][earlier]-=dot(coefficient,z);
            }
        }
    }
    reduce(finite,andOp<bool>());
    if (!finite) return false;

    if (check)
    {
        // Independent products from original coefficients, not the scatter
        // work arrays or inverse: (D+L)y-r and (D+U)z-Dy.
        const auto absoluteProduct = [](const tensor& a,const vector& x)
        {
            vector value(Zero);
            for (direction i=0;i<3;++i)
                for (direction j=0;j<3;++j)
                    value.component(i)+=mag(a.component(3*i+j)*x.component(j));
            return value;
        };
        scalar forwardError=0,backwardError=0;
        bool checkFinite=true;
        forAll(values,K)
        {
            vectorField dy(cells,Zero),dz(cells,Zero);
            vectorField scaleForward(cells,Zero),scaleBackward(cells,Zero);
            forAll(values,J)
                if (eqSystem.dVByVExists(K,J))
                {
                    const tensorField& diagonal=eqSystem.dVByV(K,J).diag();
                    for (label celli=0;celli<cells;++celli)
                    {
                        dy[celli]+=dot(diagonal[celli],forwardSolution[J][celli]);
                        dz[celli]+=dot(diagonal[celli],values[J][celli]);
                        scaleForward[celli]+=absoluteProduct(diagonal[celli],forwardSolution[J][celli]);
                        scaleBackward[celli]+=absoluteProduct(diagonal[celli],values[J][celli]);
                    }
                }
            vectorField actionForward(dy),actionBackward(dz);
            const blockFvMatrix<vector,tensor>& a=eqSystem.dVByV(K,K);
            if (a.hasLower() || a.hasUpper())
            {
                const labelUList& owner=a.lduAddr().lowerAddr();
                const labelUList& neighbour=a.lduAddr().upperAddr();
                forAll(owner,facei)
                {
                    const label own=owner[facei],nei=neighbour[facei];
                    // Direct row/column classification, independent of CSR.
                    if (own<nei)
                    {
                        actionForward[nei]+=dot(a.lower()[facei],forwardSolution[K][own]);
                        scaleForward[nei]+=absoluteProduct(a.lower()[facei],forwardSolution[K][own]);
                        actionBackward[own]+=dot(a.upper()[facei],values[K][nei]);
                        scaleBackward[own]+=absoluteProduct(a.upper()[facei],values[K][nei]);
                    }
                    else
                    {
                        actionForward[own]+=dot(a.upper()[facei],forwardSolution[K][nei]);
                        scaleForward[own]+=absoluteProduct(a.upper()[facei],forwardSolution[K][nei]);
                        actionBackward[nei]+=dot(a.lower()[facei],values[K][own]);
                        scaleBackward[nei]+=absoluteProduct(a.lower()[facei],values[K][own]);
                    }
                }
            }
            for (label celli=0;celli<cells;++celli)
                for (direction cmpt=0;cmpt<3;++cmpt)
                {
                    const scalar f=mag(actionForward[celli].component(cmpt)
                        -originalRhs[K][celli].component(cmpt))
                        /max(scaleForward[celli].component(cmpt)
                            +mag(originalRhs[K][celli].component(cmpt)),scalar(VSMALL));
                    const scalar b=mag(actionBackward[celli].component(cmpt)-dy[celli].component(cmpt))
                        /max(scaleBackward[celli].component(cmpt)+mag(dy[celli].component(cmpt)),scalar(VSMALL));
                    checkFinite=checkFinite && std::isfinite(f) && std::isfinite(b);
                    forwardError=max(forwardError,f);
                    backwardError=max(backwardError,b);
                }
        }
        reduce(checkFinite,andOp<bool>());
        reduce(forwardError,maxOp<scalar>());
        reduce(backwardError,maxOp<scalar>());
        ++nkSchurSweepChecks;
        Info<< "IHBStageSchurSweepCheck action " << nkSchurSweepChecks
            << " forwardBackwardError " << forwardError
            << " reverseBackwardError " << backwardError
            << " finite " << checkFinite << endl;
        if (!checkFinite || forwardError>1e-10 || backwardError>1e-10)
            FatalErrorInFunction << "Schur block-LUSGS triangular residual check failed"
                << exit(FatalError);
    }
    return true;
};

const auto nkApplyPressureSchur = [&]()
{
    // zU=B*rU; pressure RHS=rC-C*zU. C/G use original assembled coefficient
    // blocks and native AMI/processor actions, never physical boundary sources.
    if (!nkSchurApplyVelocityInverse(nkVectorWork)) return false;
    forAll(subTimeMeshes,K)
    {
        fvMesh& mesh=subTimeMeshes[K];
        const HBPhaseTime phaseClock
        (
            runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,
            mesh,meshPhiHB[K]
        );
        nkSchurVelocity[K].primitiveFieldRef()=nkVectorWork[K];
        nkSchurVelocity[K].correctBoundaryConditions();
        scalarField continuityAction(mesh.nCells(),Zero);
        eqSystem.dSByV(K,K).Amul(continuityAction,nkSchurVelocity[K],mesh);
        scalarField& pressureSource=nkSchurPressureMatrix[K].source();
        pressureSource=nkScalarWork[K]-continuityAction;
        bool finite=true;
        scalar maximum=0;
        forAll(pressureSource,celli)
        {
            finite=finite && std::isfinite(pressureSource[celli]);
            maximum=max(maximum,mag(pressureSource[celli]));
        }
        reduce(finite,andOp<bool>());
        reduce(maximum,maxOp<scalar>());
        if (!finite) return false;
        // Each GAMG application starts from zero, with the configured V cycles.
        // A global zero RHS is handled without a coarse-grid 0/0 normalisation.
        nkSchurPressure[K].primitiveFieldRef()=Zero;
        nkSchurPressure[K].correctBoundaryConditions();
        if (maximum>VSMALL) nkSchurPressureSolver[K].solve(nkSchurControls);
        nkScalarWork[K]=nkSchurPressure[K].primitiveField();
        forAll(nkScalarWork[K],celli)
            finite=finite && std::isfinite(nkScalarWork[K][celli]);
        reduce(finite,andOp<bool>());
        if (!finite) return false;
        eqSystem.dVByS(K,K).Amul(nkSchurGradient[K],nkSchurPressure[K],mesh);
    }
    // zU <- zU-B*G*zP, including all local cross-phase HB response in B.
    if (!nkSchurApplyVelocityInverse(nkSchurGradient)) return false;
    forAll(nkVectorWork,K) nkVectorWork[K]-=nkSchurGradient[K];
    return true;
};
