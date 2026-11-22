// nHBTurbulenceCorrectors in the caller wraps either complete update.
if (hbTurbulenceMethod=="phaseImplicit")
{
    #include "correctPhaseImplicitTurbulence.H"
}
else
{
    // Three synchronised stages give the explicit HB operator a finite imaginary
    // stability interval. A single forward-Euler/Jacobi phase sweep would amplify
    // an undamped harmonic. Spatial diffusion/destruction remain native implicit.
    if (hbSST)
        advanceHBTurbulenceStages(runTime,subTimeMeshes,snapshots,hbMotionDeltaT,
            meshPhiHB,HB,hbSSTTurbulence,turbulenceHBDamping);
    else
        advanceHBTurbulenceStages(runTime,subTimeMeshes,snapshots,hbMotionDeltaT,
            meshPhiHB,HB,hbTurbulence,turbulenceHBDamping);
}
