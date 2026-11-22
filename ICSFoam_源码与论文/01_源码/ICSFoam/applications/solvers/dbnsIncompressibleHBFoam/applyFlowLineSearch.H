// One global alpha for every phase and for p, U and the absolute face flux.
// rAU is frozen at assembly; the merit measures the current response separately.
PtrList<surfaceScalarField> flowCandidateAbsolute(nO);
scalar pressureStepMax=0,velocityStepMax=0,faceStepMax=0;
forAll(subTimeMeshes,K)
{
    const fvMesh& mesh=subTimeMeshes[K];
    const HBPhaseTime phaseClock(runTime,snapshots[K]+hbMotionDeltaT,
        hbMotionDeltaT,subTimeMeshes[K],meshPhiHB[K]);
    forAll(scalarSolIncrement[K],celli)
    {
        if (!std::isfinite(scalarSolIncrement[K][celli]))
            FatalErrorInFunction << "Nonfinite pressure increment" << exit(FatalError);
        for (direction cmpt=0;cmpt<vector::nComponents;++cmpt)
            if (!std::isfinite(vectorSolIncrement[K][celli][cmpt]))
                FatalErrorInFunction << "Nonfinite velocity increment" << exit(FatalError);
    }
    pressureStepMax=max(pressureStepMax,gMax(mag(scalarSolIncrement[K].primitiveField())));
    velocityStepMax=max(velocityStepMax,gMax(mag(vectorSolIncrement[K].primitiveField())));
    p[K] == scalarVarsPrevIter[K];
    U[K] == vectorVarsPrevIter[K];
    p[K].ref() += scalarSolIncrement[K].internalField();
    U[K].ref() += vectorSolIncrement[K].internalField();
    if (jointFlowOperator)
    {
        flowCandidateAbsolute.set(K,new surfaceScalarField
            ("qJointLineCandidate",jointFlowCandidateAbsolute[K]));
    }
    else
    {
    // Use the same candidate construction as the nonlinear stage. A second
    // pass lets inletOutlet read this candidate's flux instead of old trial phi.
    for (label boundaryPass=0;boundaryPass<2;++boundaryPass)
    {
        p[K].correctBoundaryConditions();
        U[K].correctBoundaryConditions();
        MRF[K].correctBoundaryVelocity(U[K]);
        fvScalarMatrix candidatePressureBlock
        (
            -fvm::laplacian(fvc::interpolate(flowFrozenResponse[K]),p[K])
        );
        surfaceScalarField candidateAbsolute
        (
            "qCandidatePass",fvc::flux(U[K])
          + fvc::flux(flowFrozenResponse[K]*fvc::grad(p[K]))
          + candidatePressureBlock.flux()
        );
        if (flowFacePseudoTime)
            candidateAbsolute = flowAbsolutePrevIter[K]
                +flowFaceWeight[K]*(candidateAbsolute-flowAbsolutePrevIter[K]);
        correctAbsoluteFluxBoundary(candidateAbsolute,U[K]);
        phiAbsolute[K] == candidateAbsolute;
        phi[K] == candidateAbsolute-meshPhiHB[K];
        if (boundaryPass==0)
            flowCandidateAbsolute.set
                (K,new surfaceScalarField("qLineCandidate",candidateAbsolute));
        else
            flowCandidateAbsolute[K] == candidateAbsolute;
    }
    } // legacy frozen-response candidate
    const surfaceScalarField faceStep
    (
        "lineFaceStep",mag(flowCandidateAbsolute[K]-flowAbsolutePrevIter[K])/mesh.magSf()
    );
    faceStepMax=max(faceStepMax,max(faceStep).value());
}

// Restore the COMPLETE baseline. alpha=0 must reproduce the same measurement.
forAll(subTimeMeshes,K)
{
    p[K] == scalarVarsPrevIter[K];
    U[K] == vectorVarsPrevIter[K];
    phiAbsolute[K] == flowAbsolutePrevIter[K];
    phi[K] == flowRelativePrevIter[K];
}
{
    #include "evaluateFlowResiduals.H"
}
const scalar baselineDifference=max
(
    max(mag(flowMerit-flowStartMerit)/max(1.0,flowStartMerit),
        mag(momentumMax-flowStartMomentum)/max(1.0,flowStartMomentum)),
    max(mag(continuityMax-flowStartContinuity)/max(1.0,flowStartContinuity),
        mag(fluxConsistencyMax-flowStartFlux)/max(1.0,flowStartFlux))
);
Info<< "IHBFlowAlphaZero difference " << baselineDifference
    << " merit " << flowMerit << endl;
if (!flowStateFinite || baselineDifference>1e-11)
    FatalErrorInFunction << "alpha=0 failed to restore the flow baseline"
        << exit(FatalError);

flowStepAccepted=false;
scalar acceptedAlpha=0;
scalar alpha=1.0;
for (label trial=0;trial<=flowMaxBacktracks;++trial,alpha*=0.5)
{
    forAll(subTimeMeshes,K)
    {
        const HBPhaseTime phaseClock(runTime,snapshots[K]+hbMotionDeltaT,
            hbMotionDeltaT,subTimeMeshes[K],meshPhiHB[K]);
        p[K] == scalarVarsPrevIter[K];
        U[K] == vectorVarsPrevIter[K];
        p[K].ref() += alpha*scalarSolIncrement[K].internalField();
        U[K].ref() += alpha*vectorSolIncrement[K].internalField();
        p[K].correctBoundaryConditions();
        U[K].correctBoundaryConditions();
        MRF[K].correctBoundaryVelocity(U[K]);
        phiAbsolute[K] == flowAbsolutePrevIter[K]
            + alpha*(flowCandidateAbsolute[K]-flowAbsolutePrevIter[K]);
        correctAbsoluteFluxBoundary(phiAbsolute[K],U[K]);
        phi[K] == phiAbsolute[K]-meshPhiHB[K];
    }
    {
        #include "evaluateFlowResiduals.H"
    }
    const scalar target=(1.0-flowArmijo*alpha)*flowStartMerit;
    Info<< "IHBFlowTrial alpha " << alpha << " merit " << flowMerit
        << " target " << target << " momentum " << momentumMax
        << " continuity " << continuityMax << " flux " << fluxConsistencyMax
        << " pressureStep " << alpha*pressureStepMax
        << " velocityStep " << alpha*velocityStepMax
        << " faceStep " << alpha*faceStepMax << endl;
    if (flowStateFinite && flowMerit<=target)
    {
        flowStepAccepted=true;
        acceptedAlpha=alpha;
        break;
    }
}
if (!flowStepAccepted)
{
    forAll(subTimeMeshes,K)
    {
        p[K] == scalarVarsPrevIter[K];
        U[K] == vectorVarsPrevIter[K];
        phiAbsolute[K] == flowAbsolutePrevIter[K];
        phi[K] == flowRelativePrevIter[K];
        pseudoCoField[K]=max(flowRejectFactor*pseudoCoField[K],pseudoCoNumMin);
        pseudoCoField[K].correctBoundaryConditions();
    }
    // The rejected linear residual history must not grow Co on the next retry.
    initRes.clear();
    prevRes.clear();
    {
        #include "evaluateFlowResiduals.H"
    }
    Info<< "IHBFlowRejected restoredMerit " << flowMerit
        << " pseudoCoReduction " << flowRejectFactor << endl;
}
else
{
    forAll(subTimeMeshes,K)
    {
        const HBPhaseTime phaseClock(runTime,snapshots[K]+hbMotionDeltaT,
            hbMotionDeltaT,subTimeMeshes[K],meshPhiHB[K]);
        laminarTransport[K].correct();
        if (!hbRAS) turbulence[K].correct();
    }
    Info<< "IHBFlowAccepted alpha " << acceptedAlpha
        << " merit " << flowMerit << " baselineMerit " << flowStartMerit << endl;
}
