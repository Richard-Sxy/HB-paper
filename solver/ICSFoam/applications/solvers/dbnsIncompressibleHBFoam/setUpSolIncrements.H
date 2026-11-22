// Global indexing: phase i owns scalar increment i and vector increment i.
// p is already kinematic pressure: increments are delta(p) and delta(U).
PtrList<volScalarField> scalarSolIncrement(nO);
PtrList<volVectorField> vectorSolIncrement(nO);

forAll(vectorSolIncrement, i)
{
    scalarSolIncrement.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "pIncr",
                runTime.timeName(),
                subTimeMeshes[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            subTimeMeshes[i],
            dimensionedScalar("pIncr", sqr(dimVelocity), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    vectorSolIncrement.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "UIncr",
                runTime.timeName(),
                subTimeMeshes[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            subTimeMeshes[i],
            dimensionedVector("UIncr", dimVelocity, Zero),
            zeroGradientFvPatchVectorField::typeName
        )
    );
}
