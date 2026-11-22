/* SPDX-License-Identifier: GPL-3.0-or-later
 * Independent research reference: all-phase incompressible HB JFNK.
 * OpenFOAM OpenCFD v2112. No existing user solver was used.
 */
#include "fvCFD.H"
#include "cyclicAMIPolyPatch.H"
#include "HBAlgebra.H"
#include <cmath>

using namespace Foam;

// One integrated flux per AMI source face, conservative transpose scatter.
// The geometry/weights are prescribed, never changed during a Newton solve.
void shareAMI(surfaceScalarField& flux)
{
    const polyBoundaryMesh& patches = flux.mesh().boundaryMesh();
    forAll(patches, patchi)
    {
        if (patches[patchi].type() != "cyclicAMI") continue;
        const cyclicAMIPolyPatch& pp = refCast<const cyclicAMIPolyPatch>(patches[patchi]);
        if (!pp.owner()) continue;
        const auto& address = pp.AMI().srcAddress();
        const auto& weights = pp.AMI().srcWeights();
        auto& target = flux.boundaryFieldRef()[pp.neighbPatchID()];
        target = scalar(0);
        const auto& source = flux.boundaryField()[patchi];
        forAll(address, i)
            forAll(address[i], k)
                target[address[i][k]] -= weights[i][k]*source[i];
    }
}

void shareAMI(surfaceVectorField& flux)
{
    const polyBoundaryMesh& patches = flux.mesh().boundaryMesh();
    forAll(patches, patchi)
    {
        if (patches[patchi].type() != "cyclicAMI") continue;
        const cyclicAMIPolyPatch& pp = refCast<const cyclicAMIPolyPatch>(patches[patchi]);
        if (!pp.owner()) continue;
        const auto& address = pp.AMI().srcAddress();
        const auto& weights = pp.AMI().srcWeights();
        auto& target = flux.boundaryFieldRef()[pp.neighbPatchID()];
        target = vector::zero;
        const auto& source = flux.boundaryField()[patchi];
        forAll(address, i)
        {
            vector q = source[i];
            pp.reverseTransformDirection(q, i);
            forAll(address[i], k)
                target[address[i][k]] -= weights[i][k]*q;
        }
    }
}

void checkAMI(const fvMesh& mesh, const scalar tolerance)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    forAll(patches, patchi)
    {
        if (patches[patchi].type() == "cyclicACMI"
         || patches[patchi].type() == "cyclicRepeatAMI")
            FatalErrorInFunction << "Only full-coverage cyclicAMI is supported" << exit(FatalError);
        if (patches[patchi].type() != "cyclicAMI") continue;
        const cyclicAMIPolyPatch& pp = refCast<const cyclicAMIPolyPatch>(patches[patchi]);
        if (!pp.owner()) continue;
        if (pp.applyLowWeightCorrection() || pp.createAMIFaces())
            FatalErrorInFunction << "AMI lowWeightCorrection/topology change is unsupported" << exit(FatalError);
        const auto& ami = pp.AMI();
        const auto& address = ami.srcAddress();
        const auto& weights = ami.srcWeights();
        const scalarField& area = mesh.magSf().boundaryField()[patchi];
        const scalarField& targetArea = mesh.magSf().boundaryField()[pp.neighbPatchID()];
        scalarField mappedArea(targetArea.size(), 0);
        vectorField mappedSf(targetArea.size(), vector::zero);
        scalar err = 0;
        forAll(address, i)
        {
            scalar rowSum = 0;
            forAll(address[i], k)
            {
                const label j = address[i][k];
                if (j < 0 || j >= targetArea.size() || weights[i][k] < 0)
                    FatalErrorInFunction << "Invalid serial AMI addressing/weight" << exit(FatalError);
                rowSum += weights[i][k];
                mappedArea[j] += weights[i][k]*area[i];
                vector sf = mesh.Sf().boundaryField()[patchi][i];
                pp.reverseTransformDirection(sf, i);
                mappedSf[j] -= weights[i][k]*sf;
            }
            err = max(err, mag(rowSum-1));
            err = max(err, mag(ami.srcWeightsSum()[i]-1));
        }
        forAll(targetArea, j)
        {
            err = max(err, mag(mappedArea[j]-targetArea[j])/max(targetArea[j], VSMALL));
            err = max(err, mag(ami.tgtWeightsSum()[j]-1));
            err = max(err, mag(mappedSf[j]-mesh.Sf().boundaryField()[pp.neighbPatchID()][j])
                /max(targetArea[j], VSMALL));
        }
        Info<< mesh.name() << " AMI " << pp.name() << " coverage/area error " << err << nl;
        if (err > tolerance)
            FatalErrorInFunction << "AMI coverage/area duality check failed: " << err << exit(FatalError);
    }
}

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"

    IOdictionary properties(IOobject("HBProperties", runTime.constant(), runTime,
        IOobject::MUST_READ, IOobject::NO_WRITE));
    const label harmonics = properties.get<label>("harmonics");
    const scalar omega = properties.get<scalar>("omega");
    const scalar nu = properties.get<scalar>("nu");
    const scalar uScale = properties.get<scalar>("velocityScale");
    const scalar lScale = properties.get<scalar>("lengthScale");
    const scalar rcTime = properties.get<scalar>("rcTime");
    if (harmonics < 0 || harmonics > 64 || omega <= 0 || nu <= 0 || uScale <= 0 || lScale <= 0 || rcTime <= 0
     || !std::isfinite(omega) || !std::isfinite(nu) || !std::isfinite(uScale)
     || !std::isfinite(lScale) || !std::isfinite(rcTime))
        FatalErrorInFunction << "Invalid HBProperties: positive scales, nu, omega, rcTime required" << exit(FatalError);
    const label stages = 2*harmonics+1;
    const auto D = ihb::differentiation(harmonics, omega);
    const scalar pScale = sqr(uScale);
    const scalar momentumScale = sqr(uScale)/lScale;
    const scalar continuityScale = uScale/lScale;
    const scalar amiTolerance = properties.lookupOrDefault<scalar>("amiTolerance", 1e-6);
    const scalar gclTolerance = properties.lookupOrDefault<scalar>("gclTolerance", 1e-9);
    const scalar finalTolerance = properties.lookupOrDefault<scalar>("physicalTolerance", 1e-7);
    const label refCell = properties.lookupOrDefault<label>("pRefCell", 0);
    const scalar refValue = properties.lookupOrDefault<scalar>("pRefValue", 0);
    if (!std::isfinite(amiTolerance) || amiTolerance <= 0
     || !std::isfinite(gclTolerance) || gclTolerance <= 0
     || !std::isfinite(finalTolerance) || finalTolerance <= 0 || !std::isfinite(refValue))
        FatalErrorInFunction << "Invalid finite tolerance/reference value" << exit(FatalError);

    PtrList<fvMesh> meshes(stages);
    PtrList<volVectorField> U(stages), bodyForce(stages);
    PtrList<volScalarField> p(stages);
    PtrList<surfaceScalarField> phi(stages), meshPhi(stages);
    List<bool> gauge(stages, false);
    label nCells = 0;
    for (label s=0; s<stages; ++s)
    {
        const word region("phase"+Foam::name(s));
        meshes.set(s, new fvMesh(IOobject(region, runTime.timeName(), runTime, IOobject::MUST_READ)));
        fvMesh& mesh = meshes[s];
        if (!s) nCells = mesh.nCells();
        if (mesh.nCells() != nCells || refCell < 0 || refCell >= nCells)
            FatalErrorInFunction << "All phases need identical cell correspondence and valid pRefCell" << exit(FatalError);
        if (s && (mesh.faceOwner() != meshes[0].faceOwner()
               || mesh.faceNeighbour() != meshes[0].faceNeighbour()))
            FatalErrorInFunction << "Phase mesh topology differs" << exit(FatalError);
        if (s && mesh.points() != meshes[0].points())
            FatalErrorInFunction << "Reference implementation requires fixed geometry; "
                << "rotating/ALE phase transport is specified in docs, not implemented" << exit(FatalError);
        if (s && (mesh.faces() != meshes[0].faces()
               || mesh.boundary().size() != meshes[0].boundary().size()))
            FatalErrorInFunction << "Phase face correspondence differs" << exit(FatalError);
        if (s) forAll(mesh.boundary(), patchi)
        {
            const auto& a = mesh.boundary()[patchi];
            const auto& b = meshes[0].boundary()[patchi];
            if (a.name() != b.name() || a.type() != b.type()
             || a.start() != b.start() || a.size() != b.size())
                FatalErrorInFunction << "Phase boundary correspondence differs" << exit(FatalError);
        }
        checkAMI(mesh, amiTolerance);
        U.set(s, new volVectorField(IOobject("U", runTime.timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh));
        p.set(s, new volScalarField(IOobject("p", runTime.timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh));
        if (U[s].dimensions() != dimVelocity || p[s].dimensions() != sqr(dimVelocity))
            FatalErrorInFunction << "U must be velocity and p kinematic pressure p/rho" << exit(FatalError);
        // Reject history/phi-dependent BCs: every Jv must evaluate a pure R(X).
        const wordList allowedBCs({"fixedValue", "zeroGradient", "empty", "cyclic", "cyclicAMI", "symmetryPlane", "symmetry"});
        forAll(mesh.boundary(), patchi)
        {
            if (!allowedBCs.found(U[s].boundaryField()[patchi].type())
             || !allowedBCs.found(p[s].boundaryField()[patchi].type()))
                FatalErrorInFunction << "Unsupported stateful/custom boundary on "
                    << mesh.boundary()[patchi].name() << exit(FatalError);
        }
        gauge[s] = p[s].needReference();
        if (gauge[s]) forAll(mesh.boundary(), patchi)
        {
            if (U[s].boundaryField()[patchi].type() == "zeroGradient")
                FatalErrorInFunction << "An open zeroGradient velocity boundary requires a pressure datum boundary; "
                    << "all-Neumann pressure is restricted to closed/prescribed-flux domains" << exit(FatalError);
        }
        bodyForce.set(s, new volVectorField(IOobject("bodyForce", runTime.timeName(), mesh,
            IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), mesh,
            dimensionedVector("zero", dimVelocity/dimTime, vector::zero)));
        if (bodyForce[s].dimensions() != dimVelocity/dimTime)
            FatalErrorInFunction << "bodyForce must be acceleration" << exit(FatalError);
        U[s].correctBoundaryConditions();
        p[s].correctBoundaryConditions();
        phi.set(s, new surfaceScalarField("phi", fvc::flux(U[s])));
        IOobject meshFluxIO("meshPhi", runTime.timeName(), mesh,
            IOobject::READ_IF_PRESENT, IOobject::NO_WRITE);
        const bool meshFluxSupplied = meshFluxIO.typeHeaderOk<surfaceScalarField>(true);
        meshPhi.set(s, new surfaceScalarField(meshFluxIO, mesh,
            dimensionedScalar("zero", dimVolume/dimTime, 0)));
        if (meshPhi[s].dimensions() != dimVolume/dimTime)
            FatalErrorInFunction << "meshPhi must be volumetric mesh flux" << exit(FatalError);
        scalar maxMeshFlux = gMax(mag(meshPhi[s].primitiveField()));
        forAll(mesh.boundary(), patchi)
            maxMeshFlux = max(maxMeshFlux, gMax(mag(meshPhi[s].boundaryField()[patchi])));
        if (maxMeshFlux > VSMALL)
            FatalErrorInFunction << "Nonzero meshPhi is unsupported in this fixed-geometry reference" << exit(FatalError);
        if (s && !meshFluxSupplied && mesh.points() != meshes[0].points())
            FatalErrorInFunction << "Distinct phase geometry requires explicit meshPhi at every phase" << exit(FatalError);
        // Never repair meshPhi silently: verify any AMI reconciliation below.
        surfaceScalarField meshFluxCheck("meshFluxCheck", meshPhi[s]);
        shareAMI(meshFluxCheck);
        scalar meshFluxError = 0;
        forAll(mesh.boundary(), patchi)
            meshFluxError = max(meshFluxError, gMax(mag(meshFluxCheck.boundaryField()[patchi]
                - meshPhi[s].boundaryField()[patchi])));
        if (meshFluxError > gclTolerance*uScale*sqr(lScale))
            FatalErrorInFunction << "meshPhi is not conservative across AMI" << exit(FatalError);
    }

    // Discrete GCL: div(meshPhi)_s = sum_j D_sj V_j / V_s.
    for (label s=0; s<stages; ++s)
    {
        volScalarField meshDiv(fvc::div(meshPhi[s]));
        scalar err = 0;
        for (label c=0; c<nCells; ++c)
        {
            scalar dv = 0;
            for (label j=0; j<stages; ++j) dv += D[s][j]*meshes[j].V()[c];
            err = max(err, mag(dv/meshes[s].V()[c]-meshDiv[c])/continuityScale);
        }
        Info<< "phase" << s << " normalized GCL error " << err << nl;
        if (err > gclTolerance)
            FatalErrorInFunction << "Discrete HB GCL failed" << exit(FatalError);
    }

    ihb::Vec x(stages*nCells*4, 0);
    auto index = [nCells](label s, label c, label component)
        { return 4*(s*nCells+c)+component; };
    for (label s=0; s<stages; ++s)
        for (label c=0; c<nCells; ++c)
        {
            for (label k=0; k<3; ++k) x[index(s,c,k)] = U[s][c][k]/uScale;
            x[index(s,c,3)] = p[s][c]/pScale;
        }

    const dimensionedScalar viscosity("nu", dimViscosity, nu);
    const dimensionedScalar tau("rcTime", dimTime, rcTime);
    scalar maxMomentum = GREAT, maxContinuity = GREAT, maxAMIMass = GREAT, maxAMIMomentum = GREAT;
    auto residual = [&](const ihb::Vec& state)->ihb::Vec
    {
        ihb::Vec r(state.size(), 0);
        maxMomentum = 0; maxContinuity = 0; maxAMIMass = 0; maxAMIMomentum = 0;
        // Refresh every stage before any residual evaluation, including every Jv.
        for (label s=0; s<stages; ++s)
        {
            for (label c=0; c<nCells; ++c)
            {
                for (label k=0; k<3; ++k) U[s].primitiveFieldRef()[c][k] = uScale*state[index(s,c,k)];
                p[s].primitiveFieldRef()[c] = pScale*state[index(s,c,3)];
            }
            U[s].correctBoundaryConditions();
            p[s].correctBoundaryConditions();
        }
        for (label s=0; s<stages; ++s)
        {
            fvMesh& mesh = meshes[s];
            // Share pressure traction as well as transport/diffusive momentum.
            surfaceVectorField pressureFlux("pressureFlux", fvc::interpolate(p[s])*mesh.Sf());
            shareAMI(pressureFlux);
            volVectorField gradP("gradP", fvc::surfaceIntegrate(pressureFlux));
            surfaceScalarField rawPhi("rawPhi", fvc::flux(U[s]));
            phi[s] = rawPhi + tau*(fvc::flux(gradP)-fvc::snGrad(p[s])*mesh.magSf());
            forAll(mesh.boundary(), patchi)
            {
                // Impermeable/prescribed-velocity faces must retain their flux.
                if (!U[s].boundaryField()[patchi].coupled()
                 && (U[s].boundaryField()[patchi].fixesValue()
                  || !p[s].boundaryField()[patchi].fixesValue()))
                    phi[s].boundaryFieldRef()[patchi] = rawPhi.boundaryField()[patchi];
            }
            shareAMI(phi[s]);
            surfaceScalarField relativePhi("relativePhi", phi[s]-meshPhi[s]);
            surfaceVectorField momentumFlux("momentumFlux",
                relativePhi*fvc::interpolate(U[s], relativePhi, "interpolate(U)")
                + pressureFlux - viscosity*mesh.magSf()*fvc::snGrad(U[s]));
            shareAMI(momentumFlux);
            volVectorField spatial("spatial", fvc::surfaceIntegrate(momentumFlux));
            volScalarField mass("mass", fvc::div(phi[s]));
            for (label c=0; c<nCells; ++c)
            {
                vector ru = spatial[c]-bodyForce[s][c];
                for (label j=0; j<stages; ++j)
                    ru += D[s][j]*(meshes[j].V()[c]/mesh.V()[c])*U[j][c];
                for (label k=0; k<3; ++k) r[index(s,c,k)] = ru[k]/momentumScale;
                r[index(s,c,3)] = mass[c]/continuityScale;
                maxMomentum = max(maxMomentum, mag(ru)/momentumScale);
                maxContinuity = max(maxContinuity, mag(mass[c])/continuityScale);
            }
            if (gauge[s]) r[index(s,refCell,3)] = (p[s][refCell]-refValue)/pScale;
            forAll(mesh.boundaryMesh(), patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() != "cyclicAMI") continue;
                const auto& pp = refCast<const cyclicAMIPolyPatch>(mesh.boundaryMesh()[patchi]);
                if (pp.owner())
                {
                    maxAMIMass = max(maxAMIMass,
                        mag(sum(phi[s].boundaryField()[patchi])
                          + sum(phi[s].boundaryField()[pp.neighbPatchID()]))/(uScale*sqr(lScale)));
                    vector balance = sum(momentumFlux.boundaryField()[pp.neighbPatchID()]);
                    forAll(momentumFlux.boundaryField()[patchi], facei)
                    {
                        vector q = momentumFlux.boundaryField()[patchi][facei];
                        pp.reverseTransformDirection(q, facei);
                        balance += q;
                    }
                    maxAMIMomentum = max(maxAMIMomentum, mag(balance)/(sqr(uScale)*sqr(lScale)));
                }
            }
        }
        return r;
    };

    ihb::Settings settings;
    settings.absTol = properties.lookupOrDefault<scalar>("absoluteTolerance", 1e-9);
    settings.relTol = properties.lookupOrDefault<scalar>("relativeTolerance", 1e-8);
    settings.maxNewton = properties.lookupOrDefault<label>("maxNewton", 30);
    settings.gmresRestart = properties.lookupOrDefault<label>("gmresRestart", 100);
    settings.gmresMaxIterations = properties.lookupOrDefault<label>("gmresMaxIterations", 1000);
    settings.linearRelTol = properties.lookupOrDefault<scalar>("linearRelativeTolerance", 1e-3);
    settings.fdStep = properties.lookupOrDefault<scalar>("jacobianStep", 1e-6);
    auto report = [](int iteration, double norm, double step)
    { Info<< "Newton " << iteration << " residualRMS " << norm << " step " << step << endl; };
    const ihb::Result result = ihb::newton(x, residual, settings, report);
    residual(x); // Restore accepted state after Jv probes; audit unmodified continuity too.
    Info<< "HB result: " << result.reason.c_str() << nl
        << "Newton iterations " << result.iterations << " Krylov iterations " << result.linearIterations << nl
        << "maxMomentum " << maxMomentum << " maxContinuity " << maxContinuity
        << " maxAMIMass " << maxAMIMass << " maxAMIMomentum " << maxAMIMomentum << nl;
    if (!result.converged || !std::isfinite(maxMomentum) || !std::isfinite(maxContinuity)
        || !std::isfinite(maxAMIMass) || !std::isfinite(maxAMIMomentum)
        || maxMomentum > finalTolerance || maxContinuity > finalTolerance
        || maxAMIMass > amiTolerance || maxAMIMomentum > amiTolerance)
    {
        Info<< "FAILED: no converged solution written" << endl;
        return 2;
    }
    // One output time is the entire periodic solution; phase index carries time.
    runTime.setTime(1, 1);
    for (label s=0; s<stages; ++s)
    {
        U[s].write(); p[s].write(); phi[s].write();
    }
    Info<< "CONVERGED: all-phase solution written to 1/phase*/" << endl;
    return 0;
}
