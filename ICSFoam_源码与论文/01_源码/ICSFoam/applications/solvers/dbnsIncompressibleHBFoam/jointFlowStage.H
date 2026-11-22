// Experimental Newton stage for independent cell (p,U) and face q.
// F=(V div(q), spatialMomentum(p,U,q-meshPhi)+HB(U), q-Q(p,U,q)).
// Only the velocity equation receives pseudo-time T*(U-U0). The nonlinear
// target recomputes rAU from q; assembled eqSystem is ONLY a preconditioner.
{
    struct stageFloatingPointGuard
    {
        sigFpe::ignore ignore;
        ~stageFloatingPointGuard() { std::feclearexcept(FE_ALL_EXCEPT); }
    };
    const scalar stageMomentumScale=sqr(flowVelocityScale)/flowLengthScale;
    const scalar stageContinuityScale=flowVelocityScale/flowLengthScale;
    PtrList<labelList> jointFacePatch(nO),jointFaceIndex(nO);
    PtrList<scalarField> jointFaceArea(nO),jointFaceWeight(nO),jointBaseFace(nO);
    PtrList<scalarField> jointFaceIncrement(nO),jointNegativeHq(nO);
    PtrList<scalarField> stageNegativeHc(nO),jointSavedCont(nO);
    PtrList<vectorField> stageNegativeHm(nO),jointSavedMom(nO);
    scalar jointTotalArea=0;
    bool jointGeometryValid=true;
    forAll(subTimeMeshes,K)
    {
        const fvMesh& mesh=subTimeMeshes[K];
        label count=mesh.nInternalFaces();
        forAll(mesh.boundary(),patchi)
        {
            if (isA<emptyFvPatch>(mesh.boundary()[patchi])) continue;
            if (U[K].boundaryField()[patchi].assignable())
                count+=mesh.boundary()[patchi].size();
        }
        jointFacePatch.set(K,new labelList(count,-1));
        jointFaceIndex.set(K,new labelList(count,0));
        jointFaceArea.set(K,new scalarField(count,Zero));
        jointFaceWeight.set(K,new scalarField(count,Zero));
        jointBaseFace.set(K,new scalarField(count,Zero));
        jointFaceIncrement.set(K,new scalarField(count,Zero));
        jointNegativeHq.set(K,new scalarField(count,Zero));
        stageNegativeHc.set(K,new scalarField(mesh.nCells(),Zero));
        stageNegativeHm.set(K,new vectorField(mesh.nCells(),Zero));
        jointSavedCont.set(K,new scalarField(contSource[K]));
        jointSavedMom.set(K,new vectorField(momSource[K]));
        label index=0;
        forAll(mesh.magSf().primitiveField(),facei)
        {
            const scalar area=mesh.magSf()[facei];
            jointGeometryValid=jointGeometryValid && std::isfinite(area) && area>VSMALL;
            jointFaceIndex[K][index]=facei;
            jointFaceArea[K][index]=area;
            jointFaceWeight[K][index]=area;
            jointBaseFace[K][index]=flowAbsolutePrevIter[K][facei];
            jointTotalArea+=area;
            ++index;
        }
        forAll(mesh.boundary(),patchi)
        {
            if (isA<emptyFvPatch>(mesh.boundary()[patchi])) continue;
            const scalar weight=mesh.boundary()[patchi].coupled() ? 0.5 : 1.0;
            forAll(mesh.magSf().boundaryField()[patchi],facei)
            {
                const scalar area=mesh.magSf().boundaryField()[patchi][facei];
                jointGeometryValid=jointGeometryValid && std::isfinite(area) && area>VSMALL;
                // Same global face metric as evaluateFlowResiduals.H. Fixed
                // velocity faces have zero closure error and no independent DOF.
                jointTotalArea+=weight*area;
                if (!U[K].boundaryField()[patchi].assignable()) continue;
                jointFacePatch[K][index]=patchi;
                jointFaceIndex[K][index]=facei;
                jointFaceArea[K][index]=area;
                jointFaceWeight[K][index]=weight*area;
                jointBaseFace[K][index]=flowAbsolutePrevIter[K].boundaryField()[patchi][facei];
                ++index;
            }
        }
        if (index!=count) FatalErrorInFunction << "Joint face layout mismatch" << exit(FatalError);
        // Independent traces on coupled patches. Their closure is evaluated by
        // the native phase-specific interface operators, not by local copying.
        // No assumption of AMI conservation is made: all continuity rows remain.
    }
    reduce(jointGeometryValid,andOp<bool>());
    reduce(jointTotalArea,sumOp<scalar>());
    if (!jointGeometryValid || !std::isfinite(jointTotalArea) || jointTotalArea<=VSMALL)
        FatalErrorInFunction << "Invalid joint face geometry" << exit(FatalError);

    const auto jointReadFace = [&](const surfaceScalarField& field,const label K,const label i)
    {
        const label patch=jointFacePatch[K][i],face=jointFaceIndex[K][i];
        return patch<0 ? field[face] : field.boundaryField()[patch][face];
    };
    const auto jointWriteFace = [&](surfaceScalarField& field,const label K,const label i,const scalar value)
    {
        const label patch=jointFacePatch[K][i],face=jointFaceIndex[K][i];
        if (patch<0) field.primitiveFieldRef()[face]=value;
        else field.boundaryFieldRef()[patch][face]=value;
    };
    bool stageEvaluationFinite=true;
    const auto evaluateFullStage = [&](const label iteration,const bool report=true)
    {
        stageFloatingPointGuard guard;
        bool finite=true;
        forAll(subTimeMeshes,K)
        {
            p[K] == scalarVarsPrevIter[K];
            U[K] == vectorVarsPrevIter[K];
            p[K].ref() += scalarSolIncrement[K].internalField();
            U[K].ref() += vectorSolIncrement[K].internalField();
            phiAbsolute[K] == flowAbsolutePrevIter[K];
            forAll(jointFaceIncrement[K],i)
            {
                const scalar q=jointBaseFace[K][i]
                    +jointFaceArea[K][i]*flowVelocityScale*jointFaceIncrement[K][i];
                finite=finite && std::isfinite(q);
                jointWriteFace(phiAbsolute[K],K,i,q);
            }
            phi[K] == phiAbsolute[K]-meshPhiHB[K];
            forAll(p[K],celli)
            {
                finite=finite && std::isfinite(p[K][celli]);
                for (direction c=0;c<3;++c)
                    finite=finite && std::isfinite(U[K][celli].component(c));
            }
        }
        reduce(finite,andOp<bool>());
        stageEvaluationFinite=finite;
        if (!finite) return scalar(GREAT);
        for (label pass=0;pass<2;++pass)
            forAll(subTimeMeshes,K)
            {
                const HBPhaseTime clock(runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,
                    subTimeMeshes[K],meshPhiHB[K]);
                p[K].correctBoundaryConditions();
                U[K].correctBoundaryConditions();
                MRF[K].correctBoundaryVelocity(U[K]);
                correctAbsoluteFluxBoundary(phiAbsolute[K],U[K]);
                phi[K] == phiAbsolute[K]-meshPhiHB[K];
            }
        forAll(subTimeMeshes,K) momSource[K]=Zero;
        HB.addingSource<vector>(momSource,vectorVarsGlob,0);
        hbFrame.correctSource(momSource,U,subTimeMeshes[0]);
        scalar volumeSum=0,cSum=0,mSum=0,qSum=0;
        scalar cMax=0,mMax=0,qMax=0;
        forAll(subTimeMeshes,K)
        {
            const fvMesh& mesh=subTimeMeshes[K];
            const HBPhaseTime clock(runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,
                subTimeMeshes[K],meshPhiHB[K]);
            fvVectorMatrix momentum(fvm::div(phi[K],U[K])+turbulence[K].divDevReff(U[K]));
            const volScalarField diagonal(momentum.A());
            bool validDiagonal=true;
            for (const scalar a : diagonal.primitiveField())
                validDiagonal=validDiagonal && std::isfinite(a) && a>VSMALL;
            reduce(validDiagonal,andOp<bool>());
            if (!validDiagonal)
            {
                stageEvaluationFinite=false;
                return scalar(GREAT); // Reject a probe, never divide by an invalid diagonal.
            }
            volScalarField response
            (
                IOobject("currentResponse",runTime.timeName(),mesh,
                    IOobject::NO_READ,IOobject::NO_WRITE,false),
                mesh,dimensionedScalar("zero",dimTime,0),
                zeroGradientFvPatchScalarField::typeName
            );
            response=1.0/diagonal;
            response.correctBoundaryConditions();
            const volVectorField gradP(fvc::grad(p[K]));
            fvScalarMatrix pressure(-fvm::laplacian(fvc::interpolate(response),p[K]));
            surfaceScalarField target("jointTargetFlux",fvc::flux(U[K])
                +fvc::flux(response*gradP)+pressure.flux());
            correctAbsoluteFluxBoundary(target,U[K]);
            const volScalarField divQ(fvc::div(phiAbsolute[K]));
            const volVectorField spatial(momentum.H()-diagonal*U[K]);
            stageNegativeHc[K]=-mesh.V()*divQ.primitiveField();
            stageNegativeHm[K]=mesh.V()*(spatial.primitiveField()-gradP.primitiveField())
                +momSource[K]-mesh.V()*rPseudoDeltaT[K].primitiveField()
                    *vectorSolIncrement[K].primitiveField();
            if (pressureNeedsReference[K] && pRefCells[K]>=0)
                stageNegativeHc[K][pRefCells[K]]+=pressureGaugeWeight[K]
                    *(pRefValues[K]-p[K][pRefCells[K]]);
            forAll(stageNegativeHc[K],celli)
            {
                const scalar v=mesh.V()[celli];
                const scalar c=stageNegativeHc[K][celli]/(v*stageContinuityScale);
                const scalar m=mag(stageNegativeHm[K][celli])/(v*stageMomentumScale);
                finite=finite && std::isfinite(c) && std::isfinite(m);
                volumeSum+=v;cSum+=v*sqr(c);mSum+=v*sqr(m);
                cMax=max(cMax,mag(c));mMax=max(mMax,m);
            }
            forAll(jointNegativeHq[K],i)
            {
                const scalar q=(jointReadFace(target,K,i)-jointReadFace(phiAbsolute[K],K,i))
                    /(jointFaceArea[K][i]*flowVelocityScale);
                jointNegativeHq[K][i]=q;
                finite=finite && std::isfinite(q);
                qSum+=jointFaceWeight[K][i]*sqr(q);qMax=max(qMax,mag(q));
            }
        }
        reduce(finite,andOp<bool>());reduce(volumeSum,sumOp<scalar>());
        reduce(cSum,sumOp<scalar>());reduce(mSum,sumOp<scalar>());reduce(qSum,sumOp<scalar>());
        reduce(cMax,maxOp<scalar>());reduce(mMax,maxOp<scalar>());reduce(qMax,maxOp<scalar>());
        stageEvaluationFinite=finite && volumeSum>VSMALL
            && std::isfinite(cSum) && std::isfinite(mSum) && std::isfinite(qSum);
        if (!stageEvaluationFinite) return scalar(GREAT);
        const scalar merit=Foam::sqrt((cSum+mSum)/volumeSum+qSum/jointTotalArea);
        stageEvaluationFinite=std::isfinite(merit);
        if (report) Info<< "IHBJointStage iteration " << iteration << " scaledRms " << merit
            << " continuityRms " << Foam::sqrt(cSum/volumeSum)
            << " momentumRms " << Foam::sqrt(mSum/volumeSum)
            << " closureRms " << Foam::sqrt(qSum/jointTotalArea)
            << " continuityMax " << cMax << " momentumMax " << mMax << " closureMax " << qMax << endl;
        return stageEvaluationFinite ? merit : scalar(GREAT);
    };
    scalar stageCurrentMerit=evaluateFullStage(0);
    if (!stageEvaluationFinite)
        FatalErrorInFunction << "Invalid initial joint flow stage" << exit(FatalError);
    #include "jointFlowNewtonKrylov.H"
    stageCurrentMerit=evaluateFullStage(stageNewtonSteps,false);
    if (!stageEvaluationFinite)
        FatalErrorInFunction << "Invalid final joint candidate" << exit(FatalError);
    forAll(subTimeMeshes,K)
    {
        jointFlowCandidateAbsolute.set(K,new surfaceScalarField("qJointCandidate",phiAbsolute[K]));
        p[K] == scalarVarsPrevIter[K];U[K] == vectorVarsPrevIter[K];
        phiAbsolute[K] == flowAbsolutePrevIter[K];phi[K] == flowRelativePrevIter[K];
        contSource[K]=jointSavedCont[K];momSource[K]=jointSavedMom[K];
    }
    forAll(subTimeMeshes,K)
    {
        const HBPhaseTime clock(runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,
            subTimeMeshes[K],meshPhiHB[K]);
        p[K].correctBoundaryConditions();U[K].correctBoundaryConditions();
        MRF[K].correctBoundaryVelocity(U[K]);
    }
}
