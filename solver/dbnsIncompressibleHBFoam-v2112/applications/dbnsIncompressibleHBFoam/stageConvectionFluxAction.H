// The probe state is active here. Hold its transported U fixed and
// subtract convection with the Newton-base relative flux. This isolates
// convection's face-flux action while retaining native AMI and BC terms.
forAll(subTimeMeshes,K)
{
    const HBPhaseTime phaseClock
    (
        runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,
        subTimeMeshes[K],meshPhiHB[K]
    );
    const word convectionScheme("div("+phi[K].name()+","+U[K].name()+")");
    const fvVectorMatrix convection(fvm::div(phi[K],U[K],convectionScheme));
    const fvVectorMatrix frozenConvection
    (
        fvm::div(newtonBaseFlux[K],U[K],convectionScheme)
    );
    const volVectorField convectionChange
    (
        (convection.A()-frozenConvection.A())*U[K]
        -(convection.H()-frozenConvection.H())
    );
    vectorField momentumChange(convectionChange.primitiveField());
    if (useStageAdvectiveMomentum)
    {
        const volScalarField divergenceChange(fvc::div(phi[K]-newtonBaseFlux[K]));
        vectorField rowChange(U[K].primitiveField()*divergenceChange.primitiveField());
        if (gaugeReplacesContinuity && pressureNeedsReference[K] && pRefCells[K]>=0)
            rowChange[pRefCells[K]]=Zero;
        momentumChange-=rowChange;
    }
    forAll(momentumChange,celli)
        for (direction cmpt=0;cmpt<vector::nComponents;++cmpt)
            (*convectiveFluxAction)[nkOffsets[K]+4*celli+cmpt+1]
                =momentumChange[celli].component(cmpt)/(h*stageMomentumScale);
}
