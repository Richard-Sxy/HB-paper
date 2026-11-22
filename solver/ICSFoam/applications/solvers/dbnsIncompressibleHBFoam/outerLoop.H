forAll(subTimeMeshes,K)
{
	// Reset pseudo Co number before every time step
	if (solnControl[0].dict().lookupOrDefault<bool>("resetPseudo", false) )
	{
		pseudoCoField[K] == solnControl[0].dict().lookupOrDefault<scalar>("pseudoCoNum", 1);
	}
}

scalar pseudoCoNumMin =
	solnControl[0].dict().lookupOrDefault<scalar>("pseudoCoNumMin", 0.1);

scalar pseudoCoNumMax =
	solnControl[0].dict().lookupOrDefault<scalar>("pseudoCoNumMax", 25);

scalar pseudoCoNumMaxIncr = solnControl[0].dict().lookupOrDefault<scalar>
(
	"pseudoCoNumMaxIncreaseFactor",
	1.25
);

scalar pseudoCoNumMinDecr = solnControl[0].dict().lookupOrDefault<scalar>
(
	"pseudoCoNumMinDecreaseFactor",
	0.1
);

// Constant-density continuity has no physical HB time derivative.
// Reset arrays before adding the velocity time source -V sum_k D_jk U_k.
forAll(subTimeMeshes, K)
{
    contSource[K] = Zero;
    momSource[K] = Zero;
}
// Absolute Cartesian components: template overload matches the D*I HB blocks.
HB.addingSource<vector>(momSource, vectorVarsGlob, 0);
hbFrame.correctSource(momSource,U,subTimeMeshes[0]);

scalar coNumRatio(1); // neutral until two residual histories are available

if (initRes.valid())
{
	if (!solnControl[0].firstIter() && prevRes.valid())
	{
		residualsIO& initialRes(initRes());
		residualsIO& previousRes(prevRes());

		scalar normInitSqr(0);
		scalar normPrevSqr(0);

		forAll(initialRes.vInitRes(), J)
		{
			normInitSqr += magSqr(initialRes.sInitRes()[J])
						 + magSqr(initialRes.vInitRes()[J]);

			normPrevSqr += magSqr(previousRes.sInitRes()[J])
						 + magSqr(previousRes.vInitRes()[J]);
		}

		scalar normInit = Foam::sqrt(normInitSqr);

		scalar normPrev = Foam::sqrt(normPrevSqr);

		coNumRatio = normPrev/max(normInit, VSMALL);

		coNumRatio = max(min(coNumRatio, pseudoCoNumMaxIncr), pseudoCoNumMinDecr);
	}

	prevRes.clear();
	prevRes.set(new residualsIO(initRes()));
}

forAll(subTimeMeshes,K)
{
	const fvMesh& mesh = subTimeMeshes[K];

	if (mesh.changing())
	{
		MRF[K].update();
		restoreHBMeshFlux(subTimeMeshes[K], meshPhiHB[K]);
	}

	#include "setCoAndDeltaT.H"
}

// Establish the baseline before assembly; this call does not update turbulence.
if (flowLineSearch)
{
    #include "evaluateFlowResiduals.H"
    if (!flowStateFinite)
        FatalErrorInFunction << "Invalid state at line-search baseline" << exit(FatalError);
}
const scalar flowStartMerit=flowMerit;
const scalar flowStartMomentum=momentumMax;
const scalar flowStartContinuity=continuityMax;
const scalar flowStartFlux=fluxConsistencyMax;
PtrList<surfaceScalarField> flowAbsolutePrevIter(nO),flowRelativePrevIter(nO);
PtrList<volScalarField> flowFrozenResponse(nO);
PtrList<surfaceScalarField> flowFaceWeight(nO);
scalarList pressureGaugeWeight(nO,0.0);
if (flowLineSearch)
    forAll(subTimeMeshes,K)
    {
        flowAbsolutePrevIter.set(K,new surfaceScalarField("qLineStart",phiAbsolute[K]));
        flowRelativePrevIter.set(K,new surfaceScalarField("phiLineStart",phi[K]));
    }

// Previous fields are used by the original incremental linear solver.
// Each phase has one scalar pressure and one vector velocity.
PtrList<volScalarField> scalarVarsPrevIter(nO);
PtrList<volVectorField> vectorVarsPrevIter(nO);
forAll(subTimeMeshes, i)
{
    scalarVarsPrevIter.set
    (
        i, new volScalarField(p[i].name() + "PrevIter", p[i])
    );
    vectorVarsPrevIter.set
    (
        i, new volVectorField(U[i].name() + "PrevIter", U[i])
    );
}

coupledMatrix eqSystem(subTimeMeshes[0], nO, nO, true);
// Phase K equations and their Krylov/work fields use the phase-K AMI.
eqSystem.setEquationMeshes
    (hbPhaseMeshes(subTimeMeshes), hbPhaseMeshes(subTimeMeshes));

// Matrix rows/columns are ordered as pressure scalars, then velocity vectors.
// In each phase: [ S C ; G A ] [deltaP; deltaU] = [rContinuity; rMomentum].
forAll(subTimeMeshes, J)
{
    forAll(subTimeMeshes, K)
    {
        // Matrix-free stage Newton only needs diagonal-phase p/U blocks.
        // Cross-phase pressure blocks are identically zero; keep their sparse
        // slots empty instead of storing seven zero coefficients per cell.
        // Retain the legacy full layout for solveForIncr/Picard alternatives.
        if (J==K || !stageZeroInitialGuess)
        {
            eqSystem.dSByS(J,K).diag() = Zero;
            eqSystem.dSByV(J,K).diag() = Zero;
            eqSystem.dVByS(J,K).diag() = Zero;
        }
        eqSystem.dVByV(J,K).diag() = Zero;
    }
}

forAll(subTimeMeshes, K)
{
    const fvMesh& mesh = subTimeMeshes[K];
    const HBPhaseTime phaseClock
    (
        runTime, snapshots[K] + hbMotionDeltaT, hbMotionDeltaT,
        subTimeMeshes[K], meshPhiHB[K]
    );
    volScalarField& p_ = p[K];
    volVectorField& U_ = U[K];
    MRF[K].correctBoundaryVelocity(U_);

    const bool refreshFaceFlux=!flowLineSearch;
    #include "preparePhaseEquations.H"
    // Form the predictor before insertion can modify native matrix bookkeeping.
    scalarField equationContinuityResidual(continuityResidual);
    if (flowLineSearch)
    {
        flowFrozenResponse.set(K,new volScalarField("lineFrozenResponse",rAU));
        flowFaceWeight.set(K,new surfaceScalarField
        (
            IOobject("facePseudoWeight",runTime.timeName(),mesh,
                IOobject::NO_READ,IOobject::NO_WRITE,false),
            mesh,dimensionedScalar("one",dimless,1)
        ));
        if (flowFacePseudoTime)
        {
            // (rAU_f/dtau_f)*(q-q_old) + q-Q(p,U;rAU_frozen)=0.
            // Eliminating q gives q_old + beta*(Q-q_old). beta -> 0 as
            // dtau -> 0; at a fixed point the ORIGINAL q=Q is recovered.
            flowFaceWeight[K] = 1.0/(1.0
                +fvc::interpolate(rAU)*fvc::interpolate(rPseudoDeltaT[K]));
            Info<< "IHBFacePseudoTime phase " << K << " betaMin "
                << min(flowFaceWeight[K]).value() << " betaMax "
                << max(flowFaceWeight[K]).value() << endl;
        }
        surfaceScalarField predictorPhi
        (
            "linePredictorPhi",fvc::flux(U_)
          + fvc::flux(rAU*gradP)+pressureBlock.flux()
        );
        correctAbsoluteFluxBoundary(predictorPhi,U_);
        if (flowFacePseudoTime)
        {
            predictorPhi = flowAbsolutePrevIter[K]
                +flowFaceWeight[K]*(predictorPhi-flowAbsolutePrevIter[K]);
            correctAbsoluteFluxBoundary(predictorPhi,U_);
        }
        // Eliminate the lagged face-response defect:
        // deltaq = F(p+dp,U+dU;rAU_frozen)-q, hence the continuity
        // increment RHS is -div(F(p,U;rAU_frozen)), not -div(q).
        const volScalarField predictorDiv(fvc::div(predictorPhi));
        equationContinuityResidual =
            -mesh.V()*predictorDiv.primitiveField();
    }
    else
        equationContinuityResidual = continuityResidual;
    coupledMatrix eqSystemBlock(mesh, 1, 1, true);

    // A: spatial velocity derivative plus pseudo-time diagonal V/deltaTau.
    eqSystemBlock.dVByV(0,0).insertEquation(momentumEqn);
    eqSystemBlock.dVByV(0,0).diag() +=
        mesh.V()*rPseudoDeltaT[K].primitiveField()*tensor::I;

    // G: pressure gradient; C: absolute-velocity divergence.
    // Prescribed mesh flux is independent of deltaU/deltaP, so C/S are unchanged.
    eqSystemBlock.dVByS(0,0) = fvm::grad(p_)();
    eqSystemBlock.dSByV(0,0) = flowFacePseudoTime
        ? weightedVelocityDivergence(U_,flowFaceWeight[K])()
        : fvm::div(U_)();

    // S: compact pressure-flux Jacobian approximation.
    // The full continuity residual must retain the complete face-flux formula.
    if (flowFacePseudoTime)
    {
        fvScalarMatrix stagePressureBlock
            (-fvm::laplacian(flowFaceWeight[K]*rAUf,p_));
        eqSystemBlock.dSByS(0,0).insertEquation(stagePressureBlock);
    }
    else eqSystemBlock.dSByS(0,0).insertEquation(pressureBlock);

    // insertEquation/grad/div may supply boundary sources. The full physical
    // RHS already includes boundaries, so retain it exactly once per row.
    eqSystemBlock.dSByS(0,0).source() = equationContinuityResidual;
    eqSystemBlock.dSByV(0,0).source() = Zero;
    eqSystemBlock.dVByS(0,0).source() = Zero;
    eqSystemBlock.dVByV(0,0).source() = momentumResidual;
    if (massConsistentMomentum)
    {
        // For B*dq = U*div(dq) + sum((Uf-U)*dq), continuity supplies
        // div(dq)=-div(q). Eliminate that first term using the STORED q,
        // not the scratch predictor qhat. No physical equation is changed.
        eqSystemBlock.dVByV(0,0).source() -=
            U_.primitiveField()*continuityResidual;
    }


    // Fix the constant-pressure nullspace on ONE rank per phase.
    // The unknown is deltaP, so the reference target is pRefValue - p(current).
    // Add the gauge after restoring the physical RHS; keep the continuity row.
    if (pressureNeedsReference[K] && pRefCells[K] >= 0)
    {
        const label refCell = pRefCells[K];
        const scalar refWeight = mag(eqSystemBlock.dSByS(0,0).diag()[refCell]);
        if (!std::isfinite(refWeight) || refWeight <= VSMALL)
        {
            FatalErrorInFunction
                << "Invalid pressure reference coefficient in phase " << K
                << exit(FatalError);
        }
        pressureGaugeWeight[K]=refWeight;
        eqSystemBlock.dSByS(0,0).diag()[refCell] += refWeight;
        if (gaugeReplacesContinuity)
            eqSystemBlock.dSByS(0,0).source()[refCell] = 0;
        eqSystemBlock.dSByS(0,0).source()[refCell] +=
            refWeight*(pRefValues[K] - p_[refCell]);
    }

    // Copy all same-phase spatial blocks without density/energy indexing.
    eqSystem.dSByS(K,K) += eqSystemBlock.dSByS(0,0);
    eqSystem.dSByV(K,K) = eqSystemBlock.dSByV(0,0);
    eqSystem.dVByS(K,K) = eqSystemBlock.dVByS(0,0);
    eqSystem.dVByV(K,K) += eqSystemBlock.dVByV(0,0);
}

// Physical HB couples velocity across all phases, including the diagonal.
// No HB.addBlock is applied to pressure or pressure/velocity cross-blocks.
forAll(subTimeMeshes, J)
{
    forAll(subTimeMeshes, K)
    {
        blockFvMatrix<vector,tensor> hbVelocity(subTimeMeshes[J]);
        hbVelocity.diag() = Zero;
        HB.addBlock(hbVelocity, J, K);
        hbFrame.correctBlock(hbVelocity.diag(), J, K, subTimeMeshes[J]);
        // HB contributes only cell-local diagonal coefficients. Adding the
        // whole temporary block would also add its empty interface arrays to
        // the spatial block's AMI/processor arrays.
        eqSystem.dVByV(J,K).diag() += hbVelocity.diag();
    }
}

#include "setUpSolIncrements.H"
// When only the face relation remains, use a zero p/U step and test
// the face defect itself instead of sending a zero RHS through GMRES.
const bool initialFlowSolveSkipped=flowLineSearch && (equationsConverged || stageZeroInitialGuess);
if (initialFlowSolveSkipped)
{
    initRes.clear();
    initRes.set(new residualsIO(nO,nO));
}
else
{
    // solveForIncr sums sources over ALL blocks in each equation row.
    // Save that exact r (including its single gauge contribution), and the
    // one diagonal-block source changed below. The solver does not modify M/r.
    PtrList<scalarField> pressureFluxOriginalRhs(nO);
    PtrList<vectorField> pressureFluxOriginalMomentumRhs(nO);
    PtrList<scalarField> pressureFluxDiagonalSource(nO);
    PtrList<scalarField> pressureFluxDefect(nO);
    if (pressureFluxCorrectors>0)
        forAll(subTimeMeshes,K)
        {
            const label nCells=subTimeMeshes[K].nCells();
            pressureFluxOriginalRhs.set(K,new scalarField(nCells,Zero));
            pressureFluxOriginalMomentumRhs.set(K,new vectorField(nCells,Zero));
            pressureFluxDiagonalSource.set
                (K,new scalarField(eqSystem.dSByS(K,K).source()));
            pressureFluxDefect.set(K,new scalarField(nCells,Zero));
            forAll(subTimeMeshes,J)
            {
                pressureFluxOriginalRhs[K] += eqSystem.dSByS(K,J).source();
                pressureFluxOriginalRhs[K] += eqSystem.dSByV(K,J).source();
                pressureFluxOriginalMomentumRhs[K] += eqSystem.dVByS(K,J).source();
                pressureFluxOriginalMomentumRhs[K] += eqSystem.dVByV(K,J).source();
            }
        }

    initRes.clear();
    initRes.set
    (
        new residualsIO
        (
            eqSystem.solveForIncr
            (
                scalarVarsPrevIter, vectorVarsPrevIter,
                scalarSolIncrement, vectorSolIncrement
            )
        )
    );

    if (pressureFluxCorrectors>0)
    {
        const auto updatePressureFluxDefect = [&]()
        {
            forAll(subTimeMeshes,K)
            {
                const fvMesh& mesh=subTimeMeshes[K];
                const HBPhaseTime phaseClock
                (
                    runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,
                    subTimeMeshes[K],meshPhiHB[K]
                );
                // Keep the physical name p but do not register the scratch
                // fields. This selects the SAME grad/laplacian schemes and
                // retains the native fluxRequired(p) non-orthogonal correction.
                volScalarField pressureBase
                (
                    IOobject(p[K].name(),runTime.timeName(),mesh,
                        IOobject::NO_READ,IOobject::NO_WRITE,false),
                    scalarVarsPrevIter[K]
                );
                volScalarField pressureProbe
                (
                    IOobject(p[K].name(),runTime.timeName(),mesh,
                        IOobject::NO_READ,IOobject::NO_WRITE,false),
                    scalarVarsPrevIter[K]
                );
                pressureBase.correctBoundaryConditions();
                pressureProbe.ref() += scalarSolIncrement[K].internalField();
                pressureProbe.correctBoundaryConditions();
                fvScalarMatrix basePressureBlock
                (
                    -fvm::laplacian(fvc::interpolate(flowFrozenResponse[K]),pressureBase)
                );
                fvScalarMatrix probePressureBlock
                (
                    -fvm::laplacian(fvc::interpolate(flowFrozenResponse[K]),pressureProbe)
                );
                surfaceScalarField baseFlux
                (
                    "pressureSecantBaseFlux",fvc::flux(vectorVarsPrevIter[K])
                  + fvc::flux(flowFrozenResponse[K]*fvc::grad(pressureBase))
                  + basePressureBlock.flux()
                );
                surfaceScalarField probeFlux
                (
                    "pressureSecantProbeFlux",fvc::flux(vectorVarsPrevIter[K])
                  + fvc::flux(flowFrozenResponse[K]*fvc::grad(pressureProbe))
                  + probePressureBlock.flux()
                );
                // Fixed-velocity/empty patches contribute no pressure-flux
                // increment; apply exactly the full physical boundary rule.
                correctAbsoluteFluxBoundary(baseFlux,vectorVarsPrevIter[K]);
                correctAbsoluteFluxBoundary(probeFlux,vectorVarsPrevIter[K]);
                const volScalarField fluxDifferenceDiv
                    (fvc::div(probeFlux-baseFlux));

                scalarField compactPressureAction(mesh.nCells(),Zero);
                eqSystem.dSByS(K,K).Amul
                    (compactPressureAction,scalarSolIncrement[K],mesh);
                if (pressureNeedsReference[K] && pRefCells[K]>=0)
                {
                    const label refCell=pRefCells[K];
                    compactPressureAction[refCell] -=
                        pressureGaugeWeight[K]*scalarSolIncrement[K][refCell];
                }
                // E(dp) = V div(F(U0,p0+dp)-F(U0,p0)) - S_compact dp.
                // Includes explicit non-orthogonal flux changes. With limited
                // corrected snGrad, E is a secant defect, not a linear operator.
                pressureFluxDefect[K]=
                    mesh.V()*fluxDifferenceDiv.primitiveField()
                  - compactPressureAction;
            }
        };

        const auto reportPressureFluxSecantResidual = [&](const label corrector)
        {
            PtrList<scalarField> matrixPressure(nO);
            PtrList<vectorField> matrixMomentum(nO);
            forAll(subTimeMeshes,K)
            {
                matrixPressure.set
                    (K,new scalarField(subTimeMeshes[K].nCells(),Zero));
                matrixMomentum.set
                    (K,new vectorField(subTimeMeshes[K].nCells(),Zero));
            }
            // Independent application of the assembled compact M, including
            // its AMI/processor interfaces and gauge; never solverPerf alone.
            eqSystem.matrixMul
            (
                scalarSolIncrement,vectorSolIncrement,
                matrixPressure,matrixMomentum
            );
            scalar compactSum=0,fullSum=0,rhsSum=0,maxPerVolume=0;
            scalar momentumSum=0,momentumRhsSum=0;
            bool finite=true;
            forAll(subTimeMeshes,K)
                forAll(matrixPressure[K],celli)
                {
                    const scalar compact=matrixPressure[K][celli]
                        -pressureFluxOriginalRhs[K][celli];
                    const scalar full=compact+pressureFluxDefect[K][celli];
                    const vector momentum=matrixMomentum[K][celli]
                        -pressureFluxOriginalMomentumRhs[K][celli];
                    const scalar physical=mag(full)/subTimeMeshes[K].V()[celli];
                    finite=finite && std::isfinite(compact)
                        && std::isfinite(full) && std::isfinite(physical)
                        && std::isfinite(magSqr(momentum));
                    compactSum += sqr(compact);
                    fullSum += sqr(full);
                    rhsSum += sqr(pressureFluxOriginalRhs[K][celli]);
                    momentumSum += magSqr(momentum);
                    momentumRhsSum += magSqr(pressureFluxOriginalMomentumRhs[K][celli]);
                    maxPerVolume=max(maxPerVolume,physical);
                }
            reduce(finite,andOp<bool>());
            reduce(compactSum,sumOp<scalar>());
            reduce(fullSum,sumOp<scalar>());
            reduce(rhsSum,sumOp<scalar>());
            reduce(momentumSum,sumOp<scalar>());
            reduce(momentumRhsSum,sumOp<scalar>());
            reduce(maxPerVolume,maxOp<scalar>());
            if (!finite || !std::isfinite(fullSum) || !std::isfinite(rhsSum)
             || !std::isfinite(compactSum) || !std::isfinite(momentumSum)
             || !std::isfinite(momentumRhsSum))
                FatalErrorInFunction
                    << "Nonfinite complete pressure secant residual"
                    << exit(FatalError);
            Info<< "IHBPressureFluxSecant corrector " << corrector
                << " compactPressureRelativeL2 "
                << Foam::sqrt(compactSum/max(rhsSum,VSMALL))
                << " secantPressureL2 " << Foam::sqrt(fullSum)
                << " secantPressureRelativeL2 "
                << Foam::sqrt(fullSum/max(rhsSum,VSMALL))
                << " secantPressureMaxPerVolume " << maxPerVolume
                << " momentumRelativeL2 "
                << Foam::sqrt(momentumSum/max(momentumRhsSum,VSMALL))
                << endl;
        };

        updatePressureFluxDefect();
        reportPressureFluxSecantResidual(0);
        for (label corrector=1;corrector<=pressureFluxCorrectors;++corrector)
        {
            // Secant-defect fixed point: M*x_new = r - E(x_previous). Change one source
            // block per scalar row; all momentum sources and M stay unchanged.
            forAll(subTimeMeshes,K)
                eqSystem.dSByS(K,K).source()=
                    pressureFluxDiagonalSource[K]-pressureFluxDefect[K];
            initRes.clear();
            initRes.set
            (
                new residualsIO
                (
                    eqSystem.solveForIncr
                    (
                        scalarVarsPrevIter,vectorVarsPrevIter,
                        scalarSolIncrement,vectorSolIncrement
                    )
                )
            );
            updatePressureFluxDefect();
            reportPressureFluxSecantResidual(corrector);
        }
        // Restore original bookkeeping after the diagnostic correction sweep.
        forAll(subTimeMeshes,K)
            eqSystem.dSByS(K,K).source()=pressureFluxDiagonalSource[K];
    }
}

// A pure face-closure defect can require p/U corrections even when the stored
// continuity and momentum already pass their thresholds. Keep the stage active
// after a zero-RHS initial-solve skip too.
if (fullStageCorrectors>0)
{
    #include "fullStageDefectCorrection.H"
}

forAll(subTimeMeshes,K)
{
    #include "setUpResiduals.H"
}
if (flowLineSearch)
{
    #include "applyFlowLineSearch.H"
}
else
{

forAll(subTimeMeshes, K)
{
    const HBPhaseTime phaseClock
    (
        runTime, snapshots[K] + hbMotionDeltaT, hbMotionDeltaT,
        subTimeMeshes[K], meshPhiHB[K]
    );

    // Apply deltaP/deltaU, then update transport/model coefficients.
    // phi is lagged here until the next residualsUpdate.H evaluation.
    #include "updateFields.H"
    laminarTransport[K].correct();
    if (!hbRAS) turbulence[K].correct();
}

}

// A deliberately skipped zero-RHS solve supplies no meaningful linear history.
if (initialFlowSolveSkipped)
{
    initRes.clear();
    prevRes.clear();
}
