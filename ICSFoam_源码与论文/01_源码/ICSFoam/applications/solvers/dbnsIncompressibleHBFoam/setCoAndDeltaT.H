// Local pseudo-time from convection and viscous diffusion, not sound speed.
// This is an iteration scale; it does not change the physical HB period.
const scalar maxPseudoDeltaT = solnControl[K].dict().lookupOrDefault<scalar>
(
    "maxPseudoDeltaT", 1.0
);
if (!std::isfinite(maxPseudoDeltaT) || maxPseudoDeltaT <= 0
    || !std::isfinite(pseudoCoNumMin) || !std::isfinite(pseudoCoNumMax)
    || pseudoCoNumMin <= 0 || pseudoCoNumMax < pseudoCoNumMin)
{
    FatalErrorInFunction
        << "Require maxPseudoDeltaT > 0 and 0 < pseudoCoNumMin <= pseudoCoNumMax"
        << exit(FatalError);
}

if (initRes.valid() && !solnControl[K].firstIter() && prevRes.valid())
{
    pseudoCoField[K] *= coNumRatio;
}
pseudoCoField[K] =
    max(min(pseudoCoField[K], pseudoCoNumMax), pseudoCoNumMin);
pseudoCoField[K].correctBoundaryConditions();

// nuEff is kinematic viscosity, including the selected turbulence contribution.
const volScalarField nuEff(turbulence[K].nuEff());
if (gMin(nuEff.primitiveField()) < 0)
{
    FatalErrorInFunction << "Negative effective viscosity" << exit(FatalError);
}

// Face rate has units m^3/s. Sum its positive contributions into both cells
// of an internal face; coupled boundary faces contribute to the local cell.
// phi is the current mesh-relative transport flux (lagged until evaluation).
surfaceScalarField pseudoSpectralFlux("pseudoSpectralFlux", mag(phi[K]));
pseudoSpectralFlux.setOriented(false);
surfaceScalarField viscousSpectralFlux
(
    "viscousSpectralFlux",
    2.0*fvc::interpolate(nuEff)*mesh.magSf()*mesh.deltaCoeffs()
);
viscousSpectralFlux.setOriented(false);
pseudoSpectralFlux += viscousSpectralFlux;
forAll(mesh.boundary(), patchi)
{
    if (isA<emptyFvPatch>(mesh.boundary()[patchi]))
    {
        pseudoSpectralFlux.boundaryFieldRef()[patchi] = 0.0;
    }
}

// Lambda_P = sum_f (|phi_f| + 2 nuEff_f |Sf| delta_f) / V_P.
// deltaTau_P = min(pseudoCo_P/Lambda_P, maxPseudoDeltaT).
// The inverse form also handles zero convection and zero diffusion.
const volScalarField summedRates(fvc::surfaceSum(pseudoSpectralFlux));
rPseudoDeltaT[K].ref() =
    summedRates.internalField()
   /(pseudoCoField[K].internalField()*mesh.V());
rPseudoDeltaT[K] = max
(
    rPseudoDeltaT[K],
    dimensionedScalar("minInversePseudoTime", dimless/dimTime, 1.0/maxPseudoDeltaT)
);
rPseudoDeltaT[K].correctBoundaryConditions();

Info<< "\nTime level " << K << nl
    << "Pseudo Courant No: Min: " << gMin(pseudoCoField[K].primitiveField())
    << " Mean: " << gAverage(pseudoCoField[K].primitiveField())
    << " Max: " << gMax(pseudoCoField[K].primitiveField()) << nl
    << "Pseudo deltaT: Min: " << 1.0/gMax(rPseudoDeltaT[K].primitiveField())
    << " Mean: " << gAverage(1.0/rPseudoDeltaT[K].primitiveField())
    << " Max: " << 1.0/gMin(rPseudoDeltaT[K].primitiveField()) << endl;
