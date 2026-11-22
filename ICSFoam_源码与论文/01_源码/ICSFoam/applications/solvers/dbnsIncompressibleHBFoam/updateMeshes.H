// Set all the meshes to the correct position to
// compute the grid fluxes

scalar timevalue_temp = runTime.value();
label index_temp = runTime.timeIndex();
dimensionedScalar deltaT_temp = runTime.deltaT();
const scalar hbMotionDeltaT = deltaPeriod/5000;

forAll(subTimeMeshes,K)
{
	  runTime.setTime(snapshots[K] + hbMotionDeltaT, index_temp + 1);
	  runTime.setDeltaT(hbMotionDeltaT, false);

	  subTimeMeshes[K].update();

      // Capture while the physical motion time index is still current.
      // Reading mesh.phi() after restoring runTime would clear this flux in v2112.
      if (subTimeMeshes[K].moving())
          meshPhiHB[K] = subTimeMeshes[K].phi();
      else
          meshPhiHB[K] = dimensionedScalar(dimVolume/dimTime, Zero);
      forAll(subTimeMeshes[K].boundary(), patchi)
          if (isA<emptyFvPatch>(subTimeMeshes[K].boundary()[patchi]))
              meshPhiHB[K].boundaryFieldRef()[patchi] = 0.0;
}

runTime.setTime(timevalue_temp,index_temp);
runTime.setDeltaT(deltaT_temp);

// This migration uses V*D(U): cells must retain volume and correspondence
// across HB phases (rigid rotation/translation). General deforming/topological
// meshes require a separate conservative time/geometry treatment.
forAll(subTimeMeshes, K)
{
    fvMesh& mesh = subTimeMeshes[K];
    if (mesh.nCells() != subTimeMeshes[0].nCells() || mesh.topoChanging())
        FatalErrorInFunction << "HB ALE requires unchanged cell correspondence"
            << exit(FatalError);
    const scalar volumeChange = gMax
    (
        mag(mesh.V().field()-subTimeMeshes[0].V().field())
       /subTimeMeshes[0].V().field()
    );
    const scalar probeVolumeChange = mesh.moving()
        ? gMax(mag(mesh.V().field()-mesh.V0().field())/mesh.V().field()) : 0.0;
    if (!std::isfinite(volumeChange) || !std::isfinite(probeVolumeChange)
        || max(volumeChange, probeVolumeChange) > 1e-9)
        FatalErrorInFunction
            << "Current HB ALE uses constant cell volumes; phase " << K
            << " has relative volume change " << max(volumeChange, probeVolumeChange)
            << exit(FatalError);

    restoreHBMeshFlux(mesh, meshPhiHB[K]);
    MRF[K].update();
    const HBPhaseTime phaseClock
    (
        runTime, snapshots[K] + hbMotionDeltaT, hbMotionDeltaT,
        mesh, meshPhiHB[K]
    );
    if (mesh.moving())
    {
        const surfaceScalarField schemeMeshFlux(fvc::meshPhi(U[K]));
        const scalar mismatch = max(mag(schemeMeshFlux-meshPhiHB[K])).value();
        const scalar scale = max(1.0, max(mag(meshPhiHB[K])).value());
        if (!std::isfinite(mismatch) || mismatch > 1e-12*scale)
            FatalErrorInFunction
                << "ddt(U) must supply the saved HB mesh flux; use HBM in fvSchemes"
                << exit(FatalError);
    }
    p[K].correctBoundaryConditions();
    U[K].correctBoundaryConditions();
    MRF[K].correctBoundaryVelocity(U[K]);
    // Rebuild physical mesh motion, but retain a saved independent q state.
    // Boundary correction only enforces prescribed velocity/empty-patch fluxes.
    if (!readInitialAbsoluteFlux)
        phiAbsolute[K] = fvc::flux(U[K]);
    correctAbsoluteFluxBoundary(phiAbsolute[K], U[K]);
    phi[K] = phiAbsolute[K] - meshPhiHB[K];
    Info<< "IHBALE phase " << K
        << " meshFluxMax " << max(mag(meshPhiHB[K])).value()
        << " divMeshFluxMax " << max(mag(fvc::div(meshPhiHB[K]))).value()
        << " relativeVolumeChange " << max(volumeChange, probeVolumeChange) << endl;
}

// Phase geometry must match the analytic rotor-frame rotation.
hbFrame.checkPhaseGeometry(subTimeMeshes);
