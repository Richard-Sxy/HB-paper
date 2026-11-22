// This selected probe varies U alone. Both the native residual and the
// assembled blocks use the same frozen phase flux and the same HB geometry.
for (label term=0;term<(captureVelocityMismatch ? 3 : 1);++term)
{
    velocityBlockMismatch[term].setSize(nkSize);
    velocityBlockMismatch[term]=Zero;
}
forAll(subTimeMeshes,K)
{
    fvMesh& mesh=subTimeMeshes[K];
    const HBPhaseTime phaseClock
    (
        runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,mesh,meshPhiHB[K]
    );
    const word convectionScheme("div("+phi[K].name()+","+U[K].name()+")");
    const surfaceScalarField& frozenRelativeFlux=velocityProbe==VelocityProbe::production
        ? newtonBaseFlux[K] : phi[K];
    const vectorField exactConvection
    (
        velocityTerm(fvm::div(frozenRelativeFlux,U[K],convectionScheme)(),U[K])
    );
    auto& velocity=nkSchurVelocity[K];
    vectorField& increment=velocity.primitiveFieldRef();
    forAll(increment,celli)
    {
        const label i=nkOffsets[K]+4*celli;
        increment[celli]=nkVelocityScale*vector
        (
            searchDirection[i+1],searchDirection[i+2],searchDirection[i+3]
        );
    }
    vectorField assembledConvection(mesh.nCells(),Zero);
    auditConvectionBlocks[K].Amul(assembledConvection,velocity,mesh);
    forAll(increment,celli)
    {
        const vector convectionDifference=
            (exactConvection[celli]-auditBaseConvection[K][celli])/h
            -assembledConvection[celli]/mesh.V()[celli];
        const label i=nkOffsets[K]+4*celli;
        for (direction c=0;c<3;++c)
            velocityBlockMismatch[0][i+c+1]=convectionDifference[c]/stageMomentumScale;
    }
    if (!captureVelocityMismatch) continue;

    const vectorField exactStress(velocityTerm(turbulence[K].divDevReff(U[K])(),U[K]));
    const vectorField exactHB(-momSource[K]/mesh.V());
    vectorField assembledStress(mesh.nCells(),Zero);
    auditStressBlocks[K].Amul(assembledStress,velocity,mesh);
    forAll(velocity,celli)
    {
        vector assembledHB=dot(auditHBDiagonal[K][celli],velocity[celli]);
        forAll(subTimeMeshes,J)
        {
            if (J==K || !eqSystem.dVByVExists(K,J)) continue;
            const label source=nkOffsets[J]+4*celli;
            const vector phaseIncrement=nkVelocityScale*vector
            (
                searchDirection[source+1],searchDirection[source+2],searchDirection[source+3]
            );
            assembledHB+=dot(eqSystem.dVByV(K,J).diag()[celli],phaseIncrement);
        }
        const scalar volume=mesh.V()[celli];
        const vector stressDifference=(exactStress[celli]-auditBaseStress[K][celli])/h
            -assembledStress[celli]/volume;
        const vector hbDifference=(exactHB[celli]-auditBaseTemporal[K][celli])/h
            -assembledHB/volume;
        const label i=nkOffsets[K]+4*celli;
        for (direction c=0;c<3;++c)
        {
            velocityBlockMismatch[1][i+c+1]=stressDifference[c]/stageMomentumScale;
            velocityBlockMismatch[2][i+c+1]=hbDifference[c]/stageMomentumScale;
        }
    }
}
