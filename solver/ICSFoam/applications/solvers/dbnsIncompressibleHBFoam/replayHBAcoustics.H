// Use the saved material geometry and pressure, before any flow initialisation.
fvSolution replaySolution(runTime);
const label nPhases =
    replaySolution.subDict("harmonicBalance").getOrDefault<label>("instantsNumber", 3);
if (nPhases < 2)
    FatalErrorInFunction << "HB replay needs at least two phases" << exit(FatalError);

IOdictionary motion
(
    IOobject("dynamicMeshDict", runTime.constant(), runTime,
        IOobject::READ_IF_PRESENT, IOobject::NO_WRITE, false)
);
const bool moving =
    motion.getOrDefault<word>("dynamicFvMesh", "staticFvMesh") != "staticFvMesh";

PtrList<fvMesh> meshes(nPhases);
PtrList<volScalarField> pressure(nPhases);
forAll(meshes, k)
{
    // Use the same base topology and region naming as the flow solver.
    meshes.set(k, new fvMesh(IOobject(polyMesh::defaultRegion, runTime.constant(),
        runTime, IOobject::MUST_READ, IOobject::NO_WRITE)));
    fvMesh& mesh = meshes[k];
    mesh.polyMesh::rename("subTimeLevel"+Foam::name(k));

    IOobject coords("points", runTime.timeName(), polyMesh::meshSubDir, mesh,
        IOobject::MUST_READ, IOobject::NO_WRITE, false);
    // movePoints includes parallel communication. All ranks must enter together,
    // including when an optional static-mesh point file exists on only one rank.
    const bool hasPointFile = !typeFilePath<pointIOField>(coords, false).empty();
    if (returnReduce(moving || hasPointFile, orOp<bool>()))
    {
        if (!coords.typeHeaderOk<pointIOField>(true, false))
            FatalErrorInFunction << "Missing or invalid checkpoint mesh points: "
                << coords.objectPath() << exit(FatalError);
        const pointIOField saved(coords);
        if (saved.size() != mesh.nPoints())
            FatalErrorInFunction << "Checkpoint point count mismatch in " << mesh.name()
                << exit(FatalError);
        mesh.movePoints(saved);
    }
    IOobject field("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE);
    if (!field.typeHeaderOk<volScalarField>(true, false))
        FatalErrorInFunction << "Missing or invalid checkpoint pressure: "
            << field.objectPath() << exit(FatalError);
    pressure.set(k, new volScalarField(field, mesh));
}

IOHBZoneList replayHB(meshes[0]);
replayHB.setInstants();
const scalarList& phases = replayHB.selectedSnapshots();
if (phases.size() != nPhases)
    FatalErrorInFunction << "HB phase count does not match the checkpoint" << exit(FatalError);
Info<< "IHBAcousticReplay checkpoint " << runTime.timeName() << " phases " << nPhases << endl;
runHBAcoustics(runTime, replayHB, hbMotionTimeOffset(phases[1]-phases[0]), true);
