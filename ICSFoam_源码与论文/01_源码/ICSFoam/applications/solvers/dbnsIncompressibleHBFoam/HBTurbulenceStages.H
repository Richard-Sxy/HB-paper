#ifndef incompressibleHBTurbulenceStages_H
#define incompressibleHBTurbulenceStages_H

namespace Foam
{
// Shared SSPRK3 orchestration; model classes own native spatial equations and
// the matching HB source. Physical phase clocks also select saved mesh fluxes.
template<class Model>
void advanceHBTurbulenceStages
(
    Time& time, PtrList<dynamicFvMesh>& meshes,
    const scalarList& snapshots, const scalar motionDeltaT,
    PtrList<surfaceScalarField>& meshFlux, const HBZoneList& HB,
    UPtrList<Model>& models, const scalar damping
)
{
    const label n=models.size();
    PtrList<scalarField> kStart(n), secondStart(n);
    forAll(models,j)
    {
        kStart.set(j,new scalarField(models[j].kValues()));
        secondStart.set(j,new scalarField(models[j].dissipationValues()));
    }
    for (label stage=0;stage<3;++stage)
    {
        forAll(models,j)
        {
            const HBPhaseTime clock
                (time,snapshots[j]+motionDeltaT,motionDeltaT,meshes[j],meshFlux[j]);
            models[j].prepare();
        }
        forAll(models,j) models[j].freezeTime(HB,models,j,damping);
        forAll(models,j)
        {
            const HBPhaseTime clock
                (time,snapshots[j]+motionDeltaT,motionDeltaT,meshes[j],meshFlux[j]);
            models[j].correct();
        }
        if (stage>0)
            forAll(models,j)
            {
                const HBPhaseTime clock
                    (time,snapshots[j]+motionDeltaT,motionDeltaT,meshes[j],meshFlux[j]);
                models[j].blend(stage==1 ? 0.75 : 1.0/3.0,kStart[j],secondStart[j]);
            }
    }
}

template<class Model>
void evaluateHBTurbulenceResiduals
(
    Time& time, PtrList<dynamicFvMesh>& meshes,
    const scalarList& snapshots, const scalar motionDeltaT,
    PtrList<surfaceScalarField>& meshFlux, const HBZoneList& HB,
    UPtrList<Model>& models, const scalar damping,
    scalar& kMax, scalar& secondMax, const bool audit
)
{
    forAll(models,j)
    {
        const HBPhaseTime clock
            (time,snapshots[j]+motionDeltaT,motionDeltaT,meshes[j],meshFlux[j]);
        models[j].prepare();
    }
    forAll(models,j) models[j].freezeTime(HB,models,j,damping);
    forAll(models,j)
    {
        const HBPhaseTime clock
            (time,snapshots[j]+motionDeltaT,motionDeltaT,meshes[j],meshFlux[j]);
        models[j].residuals(kMax,secondMax,audit,j);
    }
}
}
#endif
