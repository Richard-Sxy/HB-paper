// Multi-phase AMI audit for the ICSFoam coupled linear solvers.
// Builds N phase meshes exactly like dbnsIncompressibleHBFoam (one dynamicFvMesh
// per HB instant, rotor at t_K=K*T/N), assembles a block system whose diagonal
// phase blocks carry that phase's AMI interfaces, plus cell-local cross-phase
// velocity coupling, then solves it with the selected coupledMatrix solver.
// The TRUE residual is evaluated with fields owned by each phase mesh, so it
// uses the correct AMI of every phase. A solver that internally applies the
// phase-0 AMI to all phases reports convergence while this residual stays large.
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "coupledMatrix.H"
#include "blockFvOperators.H"
#include "cyclicAMIPolyPatch.H"
#include <cmath>
using namespace Foam;

static wordList phasePatchTypes(const fvMesh& mesh)
{
    wordList types(mesh.boundary().size(), "zeroGradient");
    forAll(mesh.boundary(), patchi)
    {
        const word& t = mesh.boundary()[patchi].type();
        if (polyPatch::constraintType(t)) types[patchi] = t;
    }
    return types;
}

int main(int argc, char* argv[])
{
    argList::addOption("phases", "label", "number of HB phase meshes (3)");
    argList::addOption("period", "scalar", "HB period in s (1)");
    argList::addOption("solver", "word", "GMRES or smoothSolverCoupled (GMRES)");
    argList::addOption("preconditioner", "word", "Jacobi or LUSGS (Jacobi)");
    argList::addOption("shift", "scalar", "diagonal shift 1/m^2 (2e4)");
    argList::addBoolOption("noCross", "omit cross-phase coupling");
    argList::addBoolOption("solverSequence", "mesh motion as in createMeshes.H");
    argList::addOption("maxIter", "label", "solver iterations (60)");
    argList::addBoolOption("noPhaseMeshes", "do not register phase meshes");
    #include "setRootCase.H"
    #include "createTime.H"

    const label nO = args.getOrDefault<label>("phases", 3);
    const scalar period = args.getOrDefault<scalar>("period", 1.0);

    PtrList<dynamicFvMesh> meshes(nO);
    forAll(meshes, K)
    {
        meshes.set(K, dynamicFvMesh::New(args, runTime));
        meshes[K].polyMesh::rename("subTimeLevel" + Foam::name(K));
    }
    if (args.found("solverSequence"))
    {
        // Exactly the createMeshes.H sequence of dbnsIncompressibleHBFoam.
        const scalar t0 = runTime.value();
        forAll(meshes, K)
        {
            runTime.value() = K*period/nO;
            meshes[K].update();
        }
        runTime.value() = t0;
        const label index0 = runTime.timeIndex();
        forAll(meshes, K)
        {
            runTime.setTime(K*period/nO + period/(nO*5000), index0 + 1);
            runTime.setDeltaT(period/(nO*5000), false);
            meshes[K].update();
        }
    }
    else
    forAll(meshes, K)
    {
        runTime.setTime(K*period/nO, K + 1);
        runTime.setDeltaT(period/(nO*5000), false);
        meshes[K].update();
    }
    // Report the AMI state of every phase: rotor angle seen by the interface.
    forAll(meshes, K)
        forAll(meshes[K].boundary(), patchi)
            if (isA<cyclicAMIPolyPatch>(meshes[K].boundaryMesh()[patchi]))
            {
                const auto& pp = refCast<const cyclicAMIPolyPatch>
                    (meshes[K].boundaryMesh()[patchi]);
                if (!pp.owner()) continue;
                const scalarListList& w = pp.AMI().srcWeights();
                scalar checksum = 0;
                forAll(w, i) forAll(w[i], j) checksum += (i + 1)*(j + 1)*w[i][j];
                reduce(checksum, sumOp<scalar>());
                Info<< "PhaseAMIState phase " << K << " patch " << pp.name()
                    << " firstFaceCentre "
                    << (pp.size() ? pp.faceCentres()[0] : vector::zero)
                    << " weightChecksum " << checksum << endl;
            }

    PtrList<volScalarField> p(nO), pExact(nO);
    PtrList<volVectorField> U(nO), UExact(nO);
    forAll(meshes, K)
    {
        const fvMesh& mesh = meshes[K];
        p.set(K, new volScalarField
        (
            IOobject("p", runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE, false),
            mesh, dimensionedScalar(dimless, Zero), phasePatchTypes(mesh)
        ));
        U.set(K, new volVectorField
        (
            IOobject("U", runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE, false),
            mesh, dimensionedVector(dimless, Zero), phasePatchTypes(mesh)
        ));
        pExact.set(K, new volScalarField("pExact", p[K]));
        UExact.set(K, new volVectorField("UExact", U[K]));
        // Phase-dependent, non-axisymmetric exact field, so the AMI rotor
        // angle changes the neighbour values seen across the interface.
        forAll(mesh.C(), i)
        {
            const vector& x = mesh.C()[i];
            pExact[K][i] = Foam::sin(9*x.x() + 0.7*K) + Foam::cos(7*x.y());
            UExact[K][i] = vector
            (
                Foam::cos(8*x.y() - 0.3*K), Foam::sin(6*x.x()) + 0.2*K, 0
            );
        }
        pExact[K].correctBoundaryConditions();
        UExact[K].correctBoundaryConditions();
    }

    coupledMatrix eqSystem(meshes[0], nO, nO, false);
    #ifdef ICSFOAM_PHASE_MESHES
    if (!args.found("noPhaseMeshes"))
    {
        UPtrList<const fvMesh> phaseMeshes(nO);
        forAll(meshes, K) phaseMeshes.set(K, &meshes[K]);
        eqSystem.setEquationMeshes(phaseMeshes, phaseMeshes);
    }
    #endif
    forAll(meshes, J)
        forAll(meshes, K)
        {
            // solve() reads every block source: allocate all blocks, as the
            // solver's legacy (non-stage) layout does.
            eqSystem.dSByS(J,K).diag() = Zero;
            eqSystem.dSByS(J,K).source() = Zero;
            eqSystem.dSByV(J,K).diag() = Zero;
            eqSystem.dSByV(J,K).source() = Zero;
            eqSystem.dVByS(J,K).diag() = Zero;
            eqSystem.dVByS(J,K).source() = Zero;
            eqSystem.dVByV(J,K).diag() = Zero;
            eqSystem.dVByV(J,K).source() = Zero;
        }

    const dimensionedScalar shift
    (dimless/dimArea, args.getOrDefault<scalar>("shift", 2e4));
    forAll(meshes, K)
    {
        const fvMesh& mesh = meshes[K];
        coupledMatrix block(mesh, 1, 1, false);
        fvScalarMatrix ps(fvm::Sp(shift, p[K]) - fvm::laplacian(p[K]));
        fvVectorMatrix us(fvm::Sp(shift, U[K]) - fvm::laplacian(U[K]));
        block.dSByS(0,0).diag() = Zero;
        block.dSByS(0,0).source() = Zero;
        block.dSByS(0,0).insertEquation(ps);
        block.dVByV(0,0).insertEquation(us);
        eqSystem.dSByS(K,K) += block.dSByS(0,0);
        eqSystem.dVByV(K,K) += block.dVByV(0,0);
    }
    // Cell-local cross-phase velocity coupling with an HB-like skew pattern.
    forAll(meshes, J)
        forAll(meshes, K)
            if (J != K && !args.found("noCross"))
            {
                const scalar c =
                    0.3*shift.value()
                   *Foam::sin(constant::mathematical::twoPi*(J - K)/nO);
                eqSystem.dVByV(J,K).diag() +=
                    c*meshes[J].V().field()*tensor::I;
            }

    // b = M*xExact with phase-owned fields (correct AMI per phase).
    PtrList<scalarField> bs(nO), rs(nO);
    PtrList<vectorField> bv(nO), rv(nO);
    forAll(meshes, K)
    {
        bs.set(K, new scalarField(meshes[K].nCells(), Zero));
        rs.set(K, new scalarField(meshes[K].nCells(), Zero));
        bv.set(K, new vectorField(meshes[K].nCells(), Zero));
        rv.set(K, new vectorField(meshes[K].nCells(), Zero));
    }
    eqSystem.matrixMul(pExact, UExact, bs, bv);
    forAll(meshes, J)
        forAll(meshes, K)
        {
            eqSystem.dSByS(J,K).source() = Zero;
            eqSystem.dSByV(J,K).source() = Zero;
            eqSystem.dVByS(J,K).source() = Zero;
            eqSystem.dVByV(J,K).source() = Zero;
        }
    forAll(meshes, K)
    {
        eqSystem.dSByS(K,K).source() = bs[K];
        eqSystem.dVByV(K,K).source() = bv[K];
    }

    dictionary controls, solverControls;
    const word solverName = args.getOrDefault<word>("solver", "GMRES");
    solverControls.add("preconditioner",
        args.getOrDefault<word>("preconditioner", "Jacobi"));
    solverControls.add("smoother", word("Jacobi"));
    solverControls.add("nDirections", label(40));
    solverControls.add("nSweeps", label(1));
    solverControls.add("maxIter", args.getOrDefault<label>("maxIter", 60));
    solverControls.add("tolerance", scalar(1e-14));
    solverControls.add("relTol", scalar(1e-11));
    controls.add("solver", solverName);
    controls.add(solverName, solverControls);
    eqSystem.solve(controls, p, U);

    eqSystem.matrixMul(p, U, rs, rv);
    scalar r2 = 0, b2 = 0, e2 = 0, x2 = 0;
    scalarList phaseRel(nO, 0.0);
    forAll(meshes, K)
    {
        scalar pr2 = 0, pb2 = 0;
        forAll(rs[K], i)
        {
            pr2 += sqr(rs[K][i] - bs[K][i]) + magSqr(rv[K][i] - bv[K][i]);
            pb2 += sqr(bs[K][i]) + magSqr(bv[K][i]);
            e2 += sqr(p[K][i] - pExact[K][i]) + magSqr(U[K][i] - UExact[K][i]);
            x2 += sqr(pExact[K][i]) + magSqr(UExact[K][i]);
        }
        reduce(pr2, sumOp<scalar>());
        reduce(pb2, sumOp<scalar>());
        phaseRel[K] = Foam::sqrt(pr2/max(pb2, VSMALL));
        r2 += pr2; b2 += pb2;
    }
    reduce(e2, sumOp<scalar>());
    reduce(x2, sumOp<scalar>());
    const scalar rel = Foam::sqrt(r2/max(b2, VSMALL));
    const scalar err = Foam::sqrt(e2/max(x2, VSMALL));
    Info<< "PhaseGMRESAudit ranks " << Pstream::nProcs() << " phases " << nO
        << " solver " << solverName << " trueRelativeResidual " << rel
        << " relativeSolutionError " << err << " perPhaseResidual " << phaseRel
        << endl;
    return (std::isfinite(rel) && rel < 1e-8) ? 0 : 2;
}
