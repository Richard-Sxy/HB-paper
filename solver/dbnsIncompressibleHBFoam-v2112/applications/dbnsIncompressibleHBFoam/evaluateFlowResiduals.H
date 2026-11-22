// Measure the stored (p,U,phi) state; never advance the face response.
// All phases have been updated before this evaluation. Refresh HB sources.
forAll(subTimeMeshes, K)
{
    momSource[K] = Zero;
}
// Absolute Cartesian components: template overload matches the D*I HB blocks.
HB.addingSource<vector>(momSource, vectorVarsGlob, 0);
hbFrame.correctSource(momSource,U,subTimeMeshes[0]);

momentumMax = 0;
continuityMax = 0;
// Integrated continuity defect of the gauge cell = -(net boundary/AMI flux).
scalar gaugeImbalance = 0;
referenceMax = 0;
fluxConsistencyMax = 0;
bool fieldsFinite = true;
scalar flowVolume=0, flowArea=0;
scalar momentumSum2=0, continuitySum2=0, fluxSum2=0;
const scalar momentumScale=sqr(flowVelocityScale)/flowLengthScale;
const scalar continuityScale=flowVelocityScale/flowLengthScale;
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
    const bool refreshFaceFlux=false;
    #include "preparePhaseEquations.H"

    const label gaugeCell =
        (gaugeReplacesContinuity && pressureNeedsReference[K]) ? pRefCells[K] : -1;
    if (gaugeCell >= 0)
        gaugeImbalance = max(gaugeImbalance, mag(continuityResidual[gaugeCell]));
    forAll(continuityResidual, celli)
    {
        const scalar uError = mag(momentumResidual[celli])/mesh.V()[celli];
        const scalar cError = celli == gaugeCell ? scalar(0)
            : mag(continuityResidual[celli])/mesh.V()[celli];
        fieldsFinite = fieldsFinite && std::isfinite(uError)
            && std::isfinite(cError) && std::isfinite(p_[celli]);
        for (direction cmpt=0; cmpt<vector::nComponents; ++cmpt)
            fieldsFinite = fieldsFinite && std::isfinite(U_[celli][cmpt]);
        momentumMax = max(momentumMax, uError);
        continuityMax = max(continuityMax, cError);
        const scalar cellVolume=mesh.V()[celli];
        flowVolume+=cellVolume;
        momentumSum2+=cellVolume*sqr(uError/momentumScale);
        continuitySum2+=cellVolume*sqr(cError/continuityScale);
    }
    // Optional read-only localisation of physical residuals. Report each rank's
    // maximum so that diagnostics can identify cells, interfaces and terms
    // without changing the stored pressure/velocity/face-flux state.
    if (physicalControls.lookupOrDefault<bool>("flowResidualLocationAudit",false))
    {
        const volScalarField auditMeshDivergence(fvc::div(meshPhiHB[K]));
        const auto reportLocation = [&](const word& criterion, const label celli)
        {
            if (celli<0) return;
            DynamicList<word> patches;
            forAll(mesh.boundary(),patchi)
            {
                const labelUList& cells=mesh.boundary()[patchi].faceCells();
                for (const label c : cells)
                    if (c==celli)
                    {
                        patches.append(mesh.boundary()[patchi].name());
                        break;
                    }
            }
            const scalar volume=mesh.V()[celli];
            Pout<< "IHBFlowResidualLocation phase " << K
                << " rank " << Pstream::myProcNo()
                << " criterion " << criterion << " cell " << celli
                << " centre " << mesh.C()[celli] << " volume " << volume
                << " momentum " << mag(momentumResidual[celli])/volume
                << " continuity " << mag(continuityResidual[celli])/volume
                << " temporal " << -momSource[K][celli]/volume
                << " spatial " << -spatialMomentumResidual[celli]
                << " pressureGradient " << gradP[celli]
                << " divMesh " << auditMeshDivergence[celli]
                << " patches " << patches << endl;
        };
        label momentumCell=-1,continuityCell=-1;
        scalar localMomentum=-1,localContinuity=-1;
        forAll(continuityResidual,celli)
        {
            const scalar m=mag(momentumResidual[celli])/mesh.V()[celli];
            const scalar c=mag(continuityResidual[celli])/mesh.V()[celli];
            if (m>localMomentum) { localMomentum=m; momentumCell=celli; }
            if (c>localContinuity) { localContinuity=c; continuityCell=celli; }
        }
        reportLocation("momentum",momentumCell);
        reportLocation("continuity",continuityCell);
    }
    if (pressureNeedsReference[K] && pRefCells[K] >= 0)
        referenceMax = max(referenceMax, mag(p_[pRefCells[K]]-pRefValues[K]));

    // Check the lagged response against the CURRENT convection flux.
    // Keep the same absolute/relative pair used for continuity/momentum.
    const volScalarField currentA(momentumEqn.A());
    if (gMin(currentA.primitiveField()) <= 0)
        FatalErrorInFunction << "Non-positive current momentum diagonal"
            << exit(FatalError);
    volScalarField currentResponse("currentResponse", rAU);
    currentResponse = 1.0/currentA;
    currentResponse.correctBoundaryConditions();
    fvScalarMatrix currentPressureBlock
    (
        -fvm::laplacian(fvc::interpolate(currentResponse), p_)
    );
    surfaceScalarField targetPhi
    (
        "targetPhi", fvc::flux(U_)
      + fvc::flux(currentResponse*gradP) + currentPressureBlock.flux()
    );
    if (auditFlowBoundaries) auditAbsoluteFluxBoundary(targetPhi,U_,p_);
    correctAbsoluteFluxBoundary(targetPhi, U_);
    targetPhi -= meshPhiHB[K];
    forAll(phi[K].primitiveField(), facei)
    {
        const scalar error = mag(phi[K][facei]-targetPhi[facei])/mesh.magSf()[facei];
        fieldsFinite = fieldsFinite && std::isfinite(error);
        fluxConsistencyMax = max(fluxConsistencyMax, error);
        const scalar area=mesh.magSf()[facei];
        flowArea+=area;
        fluxSum2+=area*sqr(error/flowVelocityScale);
    }
    forAll(mesh.boundary(), patchi)
    {
        const scalarField& values = phi[K].boundaryField()[patchi];
        forAll(values, facei)
        {
            fieldsFinite = fieldsFinite && std::isfinite(values[facei]);
            if (!isA<emptyFvPatch>(mesh.boundary()[patchi]))
            {
                const scalar error = mag
                (
                    values[facei]-targetPhi.boundaryField()[patchi][facei]
                )/mesh.magSf().boundaryField()[patchi][facei];
                fieldsFinite = fieldsFinite && std::isfinite(error);
                fluxConsistencyMax = max(fluxConsistencyMax, error);
                // Two processor/cyclic sides represent one coupled interface.
                const scalar area=(mesh.boundary()[patchi].coupled() ? 0.5 : 1.0)
                    *mesh.magSf().boundaryField()[patchi][facei];
                flowArea+=area;
                fluxSum2+=area*sqr(error/flowVelocityScale);
            }
        }
    }
}
auditFlowBoundaries=false;
reduce(fieldsFinite, andOp<bool>());
reduce(momentumMax, maxOp<scalar>());
reduce(continuityMax, maxOp<scalar>());
reduce(gaugeImbalance, maxOp<scalar>());
if (gaugeReplacesContinuity && gaugeImbalance > 0)
    Info<< "IHBGaugeCellImbalance maxPhaseNetFlux " << gaugeImbalance
        << " (closed-domain AMI/boundary non-conservation; not a convergence criterion)"
        << endl;
reduce(referenceMax, maxOp<scalar>());
reduce(fluxConsistencyMax, maxOp<scalar>());
reduce(flowVolume,sumOp<scalar>());
reduce(flowArea,sumOp<scalar>());
reduce(momentumSum2,sumOp<scalar>());
reduce(continuitySum2,sumOp<scalar>());
reduce(fluxSum2,sumOp<scalar>());
flowStateFinite=fieldsFinite && std::isfinite(referenceMax)
    && std::isfinite(momentumSum2) && std::isfinite(continuitySum2)
    && std::isfinite(fluxSum2) && flowVolume>0 && flowArea>0;
if (!flowStateFinite && !flowLineSearch)
    FatalErrorInFunction << "Nonfinite physical fields or residuals" << exit(FatalError);
flowMomentumRms=flowStateFinite ? Foam::sqrt(momentumSum2/flowVolume) : GREAT;
flowContinuityRms=flowStateFinite ? Foam::sqrt(continuitySum2/flowVolume) : GREAT;
flowFluxRms=flowStateFinite ? Foam::sqrt(fluxSum2/flowArea) : GREAT;
flowMerit=flowStateFinite ? Foam::sqrt
(
    sqr(flowMomentumRms)+sqr(flowContinuityRms)+sqr(flowFluxRms)
  + sqr(referenceMax/sqr(flowVelocityScale))
) : GREAT;

equationsConverged = flowStateFinite && momentumMax <= momentumTolerance
    && continuityMax <= continuityTolerance && referenceMax <= referenceTolerance;
physicalConverged = equationsConverged && fluxConsistencyMax <= fluxTolerance;
Info<< "IHBPhysicalResidual momentum " << momentumMax
    << " continuity " << continuityMax << " reference " << referenceMax
    << " fluxConsistency " << fluxConsistencyMax << endl;

Info<< "IHBFlowMerit " << flowMerit
    << " momentumRms " << flowMomentumRms
    << " continuityRms " << flowContinuityRms
    << " fluxRms " << flowFluxRms << endl;
