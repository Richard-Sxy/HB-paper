// Each phase owns one scalar (kinematic pressure p) and one vector (U).
// Keep the original ICSFoam phase/variable layout for later block assembly.
PtrList<PtrList<volScalarField>> scalarVarsGlob(nO);
PtrList<PtrList<volVectorField>> vectorVarsGlob(nO);

// Non-owning views: p[i] and U[i] refer to the same fields held above.
UPtrList<volScalarField> p(nO);
UPtrList<volVectorField> U(nO);

// One pressure gauge per phase. The reference cell is local to its owner rank;
// setRefCell returns -1 on all other ranks, or if no reference is needed.
labelList pRefCells(nO, -1);
scalarList pRefValues(nO, 0.0);
boolList pressureNeedsReference(nO, false);

PtrList<scalarField> contSource(nO);
PtrList<vectorField> momSource(nO);
// phi is mesh-relative transport flux; phiAbsolute enforces div(U)=0.
PtrList<surfaceScalarField> phi(nO);
PtrList<surfaceScalarField> phiAbsolute(nO);
PtrList<surfaceScalarField> meshPhiHB(nO);

// A line-search restart must retain the independent absolute face-flux state.
// Check only the requested time, never a previous instance of phiAbsolute.
boolList phaseAbsoluteFluxRead(nO, false);
bool readInitialAbsoluteFlux=false;

// Construct transport before turbulence so it outlives the turbulence model.
PtrList<singlePhaseTransportModel> laminarTransport(nO);
PtrList<incompressible::turbulenceModel> turbulence(nO);
UPtrList<HBKEpsilon> hbTurbulence(nO);
UPtrList<HBKOmegaSST> hbSSTTurbulence(nO);
bool hbRAS=false;
bool hbSST=false;
word hbRASModel;

PtrList<volScalarField> pseudoCoField(nO);
PtrList<volScalarField> rPseudoDeltaT(nO);
autoPtr<residualsIO> prevRes(nullptr);
autoPtr<residualsIO> initRes(nullptr);
PtrList<IOMRFCoupledZoneList> MRF(nO);

bool anyInitialAbsoluteFlux=false;
bool allInitialAbsoluteFlux=true;
forAll(subTimeMeshes, i)
{
    IOobject absoluteFluxHeader
    (
        "phiAbsolute", runTime.timeName(), subTimeMeshes[i],
        IOobject::READ_IF_PRESENT, IOobject::NO_WRITE, false
    );
    // Check header existence before class: a fresh IOobject itself has a
    // non-empty default headerClassName even when no field file exists.
    phaseAbsoluteFluxRead[i] =
        absoluteFluxHeader.typeHeaderOk<surfaceScalarField>(false, false);
    if (phaseAbsoluteFluxRead[i]
        && absoluteFluxHeader.headerClassName() != surfaceScalarField::typeName)
        FatalErrorInFunction
            << "Invalid phiAbsolute field class in phase " << i
            << " at time " << runTime.timeName()
            << "; expected surfaceScalarField" << exit(FatalError);
    anyInitialAbsoluteFlux = anyInitialAbsoluteFlux || phaseAbsoluteFluxRead[i];
    allInitialAbsoluteFlux = allInitialAbsoluteFlux && phaseAbsoluteFluxRead[i];
}
reduce(anyInitialAbsoluteFlux, orOp<bool>());
reduce(allInitialAbsoluteFlux, andOp<bool>());
if (anyInitialAbsoluteFlux && !allInitialAbsoluteFlux)
    FatalErrorInFunction
        << "phiAbsolute must be present for every HB phase and every rank,"
        << " or absent everywhere, at the selected time " << runTime.timeName()
        << exit(FatalError);
readInitialAbsoluteFlux=allInitialAbsoluteFlux;
Info<< "IHBInitialAbsoluteFlux read=" << readInitialAbsoluteFlux
    << " time=" << runTime.timeName() << endl;

Info<< "Reading incompressible fields p and U" << nl << endl;
forAll(subTimeMeshes, i)
{
    scalarVarsGlob.set(i, new PtrList<volScalarField>(1));
    vectorVarsGlob.set(i, new PtrList<volVectorField>(1));

    scalarVarsGlob[i].set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                "p", runTime.timeName(), subTimeMeshes[i],
                IOobject::MUST_READ, IOobject::AUTO_WRITE
            ),
            subTimeMeshes[i]
        )
    );
    vectorVarsGlob[i].set
    (
        0,
        new volVectorField
        (
            IOobject
            (
                "U", runTime.timeName(), subTimeMeshes[i],
                IOobject::MUST_READ, IOobject::AUTO_WRITE
            ),
            subTimeMeshes[i]
        )
    );

    p.set(i, &scalarVarsGlob[i][0]);
    U.set(i, &vectorVarsGlob[i][0]);

    // p means pressure divided by the constant reference density.
    if (p[i].dimensions() != sqr(dimVelocity)
        || U[i].dimensions() != dimVelocity)
    {
        FatalErrorInFunction
            << "Expected kinematic pressure p [0 2 -2 0 0 0 0] and velocity U"
            << exit(FatalError);
    }
    subTimeMeshes[i].setFluxRequired(p[i].name());

    // Native OpenFOAM selection supports pRefCell or pRefPoint plus pRefValue.
    // It checks needReference() globally and selects only one owning rank.
    pressureNeedsReference[i] = setRefCell
    (
        p[i], solnControl[i].dict(), pRefCells[i], pRefValues[i]
    );
    if (pressureNeedsReference[i] && !std::isfinite(pRefValues[i]))
    {
        FatalErrorInFunction << "pRefValue must be finite" << exit(FatalError);
    }

    // phi is volume flux, not the original compressible mass flux.
    // Initial value only; updateMeshes.H rebuilds it after capturing mesh motion.
    phi.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi", runTime.timeName(), subTimeMeshes[i],
                IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
            ),
            fvc::flux(U[i])
        )
    );
    if (phi[i].dimensions() != dimVolume/dimTime)
    {
        FatalErrorInFunction
            << "Expected volume flux phi; a compressible mass-flux file cannot be reused"
            << exit(FatalError);
    }

    phiAbsolute.set
    (
        i, new surfaceScalarField
        (
            IOobject
            (
                "phiAbsolute", runTime.timeName(), subTimeMeshes[i],
                phaseAbsoluteFluxRead[i] ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::flux(U[i])
        )
    );
    if (phiAbsolute[i].dimensions() != dimVolume/dimTime)
        FatalErrorInFunction
            << "Expected volume flux phiAbsolute [0 3 -1 0 0 0 0] in phase " << i
            << "; a compressible mass-flux file cannot be reused"
            << exit(FatalError);

    meshPhiHB.set
    (
        i, new surfaceScalarField
        (
            IOobject("meshPhiHB", runTime.timeName(), subTimeMeshes[i],
                     IOobject::NO_READ, IOobject::AUTO_WRITE),
            subTimeMeshes[i], dimensionedScalar(dimVolume/dimTime, Zero)
        )
    );
    meshPhiHB[i].setOriented();

    contSource.set(i, new scalarField(subTimeMeshes[i].nCells(), Zero));
    momSource.set(i, new vectorField(subTimeMeshes[i].nCells(), Zero));
}

Info<< "Creating incompressible transport and turbulence models" << nl << endl;
forAll(subTimeMeshes, i)
{
    laminarTransport.set(i, new singlePhaseTransportModel(U[i], phi[i]));
    IOdictionary modelDict
    (
        IOobject("turbulenceProperties", runTime.constant(), subTimeMeshes[i],
                 IOobject::MUST_READ, IOobject::NO_WRITE, false)
    );
    const word simulationType(modelDict.get<word>("simulationType"));
    const bool useRAS=simulationType=="RAS";
    if (i==0) hbRAS=useRAS;
    if (useRAS!=hbRAS || (simulationType!="laminar" && !useRAS))
        FatalErrorInFunction << "Use laminar or supported RAS consistently in all HB phases"
            << exit(FatalError);
    if (useRAS)
    {
        const word modelName(modelDict.subDict("RAS").get<word>("RASModel"));
        if (i==0)
        {
            hbRASModel=modelName;
            hbSST=modelName=="kOmegaSST";
        }
        if (modelName!=hbRASModel)
            FatalErrorInFunction << "Use the same RAS model in every HB phase"
                << exit(FatalError);
        bool active=false;
        if (modelName=="kEpsilon")
        {
            HBKEpsilon* model=new HBKEpsilon(U[i],phi[i],laminarTransport[i]);
            turbulence.set(i,model);
            hbTurbulence.set(i,model);
            active=model->active();
        }
        else if (hbSST)
        {
            HBKOmegaSST* model=new HBKOmegaSST(U[i],phi[i],laminarTransport[i]);
            turbulence.set(i,model);
            hbSSTTurbulence.set(i,model);
            active=model->active();
        }
        else
            FatalErrorInFunction << "Supported HB RAS models: kEpsilon, kOmegaSST"
                << exit(FatalError);
        if (!active)
            FatalErrorInFunction << "Set RAS turbulence on, or select simulationType laminar"
                << exit(FatalError);
    }
    else
        turbulence.set(i,incompressible::turbulenceModel::New(U[i],phi[i],laminarTransport[i]));
}

// Retain the original pseudo-time storage for the incompressible local step.
forAll(rPseudoDeltaT, i)
{
    rPseudoDeltaT.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rPseudoDeltaT", runTime.timeName(), subTimeMeshes[i],
                IOobject::NO_READ, IOobject::NO_WRITE
            ),
            subTimeMeshes[i],
            dimensionedScalar("rPseudoDeltaT", dimless/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    pseudoCoField.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "pseudoCoField", runTime.timeName(), subTimeMeshes[i],
                IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
            ),
            subTimeMeshes[i],
            dimensionedScalar
            (
                "pseudoCoNumInit", dimless,
                solnControl[0].dict().lookupOrDefault<scalar>("pseudoCoNum", 1)
            )
        )
    );
}

forAll(MRF, i)
{
    MRF.set(i, new IOMRFCoupledZoneList(subTimeMeshes[i]));
    MRF[i].correctBoundaryVelocity(U[i]);
}
