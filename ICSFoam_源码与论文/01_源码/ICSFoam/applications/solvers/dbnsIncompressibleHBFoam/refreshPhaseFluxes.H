// Explicitly advance the lagged face-response relation once. Keeping this
// separate makes repeated physical residual measurements reproducible.
forAll(subTimeMeshes, K)
{
    const fvMesh& mesh = subTimeMeshes[K];
    volScalarField& p_ = p[K];
    volVectorField& U_ = U[K];
    const HBPhaseTime phaseClock
    (
        runTime, snapshots[K] + hbMotionDeltaT, hbMotionDeltaT,
        subTimeMeshes[K], meshPhiHB[K]
    );
    p_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    MRF[K].correctBoundaryVelocity(U_);
    const bool refreshFaceFlux=true;
    #include "preparePhaseEquations.H"
}
