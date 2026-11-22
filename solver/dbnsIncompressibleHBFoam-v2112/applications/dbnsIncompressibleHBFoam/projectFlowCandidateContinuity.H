// Apply the same pressure correction to p, U and q. In the frozen face
// response the velocity correction cancels the reconstructed cell gradient,
// leaving the native pressure-matrix face flux as the change in q.
const scalar eta=runTime.timeIndex()>=stageContinuityStartIndex
    ? stageContinuityEta : scalar(1);
const volScalarField targetDiv((1-eta)*fvc::div(flowAbsolutePrevIter[K]));
const volScalarField initialDefect
(
    targetDiv-fvc::div(flowCandidateAbsolute[K])
);
const scalar before=gMax(mag(initialDefect.primitiveField()));
if (before>continuityTolerance)
{
    wordList pressureTypes(p[K].boundaryField().types());
    for (const fvPatchScalarField& patch : p[K].boundaryField())
        if (!patch.coupled() && patch.type()!="empty"
         && patch.type()!="fixedValue" && patch.type()!="zeroGradient")
            FatalErrorInFunction << "Unsupported pressure patch for candidate projection: "
                << patch.patch().name() << " (" << patch.type() << ')'
                << exit(FatalError);

    volScalarField correction
    (
        IOobject("pCandidateProjection",runTime.timeName(),mesh,
            IOobject::NO_READ,IOobject::NO_WRITE,false),
        mesh,dimensionedScalar("zero",p[K].dimensions(),0),pressureTypes
    );
    mesh.setFluxRequired(correction.name());
    surfaceScalarField response(fvc::interpolate(flowFrozenResponse[K]));
    const word laplacianName("laplacian("+response.name()+","+p[K].name()+")");
    forAll(mesh.boundary(),patchi)
        if (!U[K].boundaryField()[patchi].assignable()
         || isA<emptyFvPatch>(mesh.boundary()[patchi]))
            response.boundaryFieldRef()[patchi]=0;

    dictionary controls;
    controls.set("solver",word("GAMG"));
    controls.set("smoother",word("GaussSeidel"));
    controls.set("agglomerator",word("faceAreaPair"));
    controls.set("directSolveCoarsest",false);
    controls.set("tolerance",scalar(1e-12));
    controls.set("relTol",scalar(1e-5));
    controls.set("maxIter",label(30));
    controls.set("log",label(0));
    surfaceScalarField projectedQ("qCandidateProjected",flowCandidateAbsolute[K]);
    for (label pass=0;pass<2;++pass)
    {
        fvScalarMatrix projection
        (
            -fvm::laplacian(response,correction,laplacianName)==initialDefect
        );
        if (pressureNeedsReference[K] && pRefCells[K]>=0)
            projection.setReference(pRefCells[K],0);
        projection.solve(controls);
        if (pass==1)
            projectedQ += projection.flux();
    }
    const volScalarField remainingDefect(targetDiv-fvc::div(projectedQ));
    const scalar after=gMax(mag(remainingDefect.primitiveField()));
    const word gradientName("grad("+p[K].name()+")");
    const volVectorField velocityCorrection
    (
        -flowFrozenResponse[K]*fvc::grad(correction,gradientName)
    );
    bool finite=std::isfinite(after);
    forAll(correction,celli)
    {
        finite=finite && std::isfinite(correction[celli]);
        for (direction cmpt=0;cmpt<vector::nComponents;++cmpt)
            finite=finite && std::isfinite(velocityCorrection[celli][cmpt]);
    }
    reduce(finite,andOp<bool>());
    const bool retain=finite && after<before;
    if (retain)
    {
        scalarSolIncrement[K].primitiveFieldRef()+=correction.primitiveField();
        vectorSolIncrement[K].primitiveFieldRef()+=velocityCorrection.primitiveField();
        flowCandidateAbsolute[K] == projectedQ;
    }
    Info<< "IHBFlowCandidateProjection phase " << K
        << " before " << before << " after " << after << " used " << retain << endl;
}
