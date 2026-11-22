#include "evaluateFlowResiduals.H"

// Native spatial terms + HB time terms; no relaxation or pseudo-time diagonal.
kResidualMax=0; epsilonResidualMax=0;
if (hbRAS && !freezeTurbulence)
{
    if (hbSST)
        evaluateHBTurbulenceResiduals(runTime,subTimeMeshes,snapshots,hbMotionDeltaT,
            meshPhiHB,HB,hbSSTTurbulence,turbulenceHBDamping,
            kResidualMax,epsilonResidualMax,turbulenceResidualAudit);
    else
        evaluateHBTurbulenceResiduals(runTime,subTimeMeshes,snapshots,hbMotionDeltaT,
            meshPhiHB,HB,hbTurbulence,turbulenceHBDamping,
            kResidualMax,epsilonResidualMax,turbulenceResidualAudit);
    reduce(kResidualMax,maxOp<scalar>());
    reduce(epsilonResidualMax,maxOp<scalar>());
    Info<< "IHBTurbulenceResidual k " << kResidualMax
        << (hbSST ? " omega " : " epsilon ") << epsilonResidualMax << endl;
}
if (hbRAS && freezeTurbulence)
{
    // A frozen-model diagnostic is never labelled a converged RANS-HB solve.
    kResidualMax=GREAT; epsilonResidualMax=GREAT;
    Info<< "IHBTurbulenceFrozen diagnosticOnly=1" << endl;
}
turbulenceConverged=kResidualMax<=kTolerance
    && epsilonResidualMax<=secondTurbulenceTolerance;
physicalConverged=physicalConverged && turbulenceConverged;
