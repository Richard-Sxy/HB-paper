// Frozen-M nonlinear stage correction. Included after the first coupled solve.
// x is ALWAYS a total increment from the same outer baseline, never accumulated
// on the previous candidate. H = (Fc + gauge, Fm + T*deltaU), T = V/deltaTau.
// Each correction solves M*xNew = M*x - H(x). All native interfaces remain in M.
// nut, rAU, T, M, the gauge weight and merit scales remain fixed here.
{
    // Clear exceptions raised by a rejected probe before restoring native traps.
    // Member destruction follows this destructor body, so the ordering matters.
    struct stageFloatingPointGuard
    {
        sigFpe::ignore ignore;
        ~stageFloatingPointGuard() { std::feclearexcept(FE_ALL_EXCEPT); }
    };
    const label nBlocks=nO*nO;
    PtrList<scalarField> stageSavedSS(nBlocks),stageSavedSV(nBlocks);
    PtrList<vectorField> stageSavedVS(nBlocks),stageSavedVV(nBlocks);
    PtrList<scalarField> stageSavedContinuitySource(nO);
    PtrList<vectorField> stageSavedMomentumSource(nO);
    PtrList<scalarField> stageNegativeHc(nO),stageMatrixPressure(nO);
    PtrList<vectorField> stageNegativeHm(nO),stageMatrixMomentum(nO);
    forAll(subTimeMeshes,K)
    {
        const label nCells=subTimeMeshes[K].nCells();
        stageSavedContinuitySource.set(K,new scalarField(contSource[K]));
        stageSavedMomentumSource.set(K,new vectorField(momSource[K]));
        stageNegativeHc.set(K,new scalarField(nCells,Zero));
        stageNegativeHm.set(K,new vectorField(nCells,Zero));
        stageMatrixPressure.set(K,new scalarField(nCells,Zero));
        stageMatrixMomentum.set(K,new vectorField(nCells,Zero));
        // Newton only applies coefficients; it never changes eqSystem sources.
        // Save the N^2 source arrays solely for the Picard alternative below.
        if (stageNewtonSteps==0)
        {
            forAll(subTimeMeshes,J)
            {
                const label block=K*nO+J;
                stageSavedSS.set(block,new scalarField(eqSystem.dSByS(K,J).source()));
                stageSavedSV.set(block,new scalarField(eqSystem.dSByV(K,J).source()));
                stageSavedVS.set(block,new vectorField(eqSystem.dVByS(K,J).source()));
                stageSavedVV.set(block,new vectorField(eqSystem.dVByV(K,J).source()));
            }
        }
    }

    const scalar stageMomentumScale=sqr(flowVelocityScale)/flowLengthScale;
    const scalar stageContinuityScale=
        (flowVelocityScale/flowLengthScale)/stageContinuityWeight;
    const scalar eta = runTime.timeIndex()>=stageContinuityStartIndex
        ? stageContinuityEta : scalar(1);
    PtrList<scalarField> stageContinuityOffset(nO);
    if (eta<1)
    {
        forAll(subTimeMeshes,K)
        {
            const volScalarField anchorDiv(fvc::div(flowAbsolutePrevIter[K]));
            stageContinuityOffset.set
            (
                K,
                new scalarField
                (
                    (1-eta)*subTimeMeshes[K].V()*anchorDiv.primitiveField()
                )
            );
        }
    }
    Info<< "IHBStageContinuity time " << runTime.timeName()
        << " eta " << eta << " normWeight " << stageContinuityWeight
        << " momentumPseudoTime " << useStageMomentumPseudoTime << endl;
    if (stageFluxRestoreFraction<1)
        Info<< "IHBStageFluxRestore time " << runTime.timeName()
            << " fraction " << stageFluxRestoreFraction << endl;
    bool stageEvaluationFinite=true;
    const auto evaluateFullStage = [&]
    (
        const label corrector,
        const bool report=true,
        const PtrList<surfaceScalarField>* fixedAbsoluteFlux=nullptr
    )
    {
        // A trial overflow is a rejected nonlinear candidate, not a process
        // failure. Restore the previous SIGFPE setting on every return path.
        stageFloatingPointGuard stageFpeGuard;
        stageEvaluationFinite=true;
        forAll(subTimeMeshes,K)
            forAll(scalarSolIncrement[K],celli)
            {
                stageEvaluationFinite=stageEvaluationFinite
                    && std::isfinite(scalarSolIncrement[K][celli]);
                for (direction cmpt=0;cmpt<vector::nComponents;++cmpt)
                    stageEvaluationFinite=stageEvaluationFinite
                        && std::isfinite(vectorSolIncrement[K][celli][cmpt]);
            }
        reduce(stageEvaluationFinite,andOp<bool>());
        if (!stageEvaluationFinite)
        {
            if (report) Info<< "IHBFullStage corrector " << corrector
                << " finite 0 scaledRms " << GREAT << endl;
            return scalar(GREAT);
        }
        // Restore the same flux before resetting flux-dependent patch branches.
        // Assign ALL phase internals before evaluating any HB temporal source.
        forAll(subTimeMeshes,K)
        {
            p[K] == scalarVarsPrevIter[K];
            U[K] == vectorVarsPrevIter[K];
            phiAbsolute[K] == flowAbsolutePrevIter[K];
            phi[K] == flowRelativePrevIter[K];
            if (fixedAbsoluteFlux)
            {
                phiAbsolute[K] == (*fixedAbsoluteFlux)[K];
                phi[K] == (*fixedAbsoluteFlux)[K]-meshPhiHB[K];
            }
            p[K].ref() += scalarSolIncrement[K].internalField();
            U[K].ref() += vectorSolIncrement[K].internalField();
        }
        // Baseline plus a finite increment can itself overflow. Check before
        // native boundary/gradient operations and make the decision collective.
        forAll(subTimeMeshes,K)
            forAll(p[K],celli)
            {
                stageEvaluationFinite=stageEvaluationFinite
                    && std::isfinite(p[K][celli]);
                for (direction cmpt=0;cmpt<vector::nComponents;++cmpt)
                    stageEvaluationFinite=stageEvaluationFinite
                        && std::isfinite(U[K][celli][cmpt]);
            }
        reduce(stageEvaluationFinite,andOp<bool>());
        if (!stageEvaluationFinite)
        {
            if (report) Info<< "IHBFullStage corrector " << corrector
                << " finite 0 scaledRms " << GREAT << endl;
            return scalar(GREAT);
        }
        // The second pass lets flux-dependent patches read the candidate phi
        // instead of inheriting the previous inner candidate's boundary branch.
        for (label boundaryPass=0;boundaryPass<2;++boundaryPass)
            forAll(subTimeMeshes,K)
            {
                const HBPhaseTime phaseClock
                (
                    runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,
                    subTimeMeshes[K],meshPhiHB[K]
                );
                p[K].correctBoundaryConditions();
                U[K].correctBoundaryConditions();
                MRF[K].correctBoundaryVelocity(U[K]);
                if (fixedAbsoluteFlux) continue;
                // Keep the physical p name: grad(p), fluxRequired(p) and the
                // native explicit non-orthogonal pressure flux must be retained.
                fvScalarMatrix stagePressureBlock
                (
                    -fvm::laplacian(fvc::interpolate(flowFrozenResponse[K]),p[K])
                );
                phiAbsolute[K] == fvc::flux(U[K])
                    + fvc::flux(flowFrozenResponse[K]*fvc::grad(p[K]))
                    + stagePressureBlock.flux();
                if (stageFluxClosureOffset.set(K))
                    phiAbsolute[K] += stageFluxClosureOffset[K];
                if (flowFacePseudoTime)
                    phiAbsolute[K] == flowAbsolutePrevIter[K]
                        +flowFaceWeight[K]*(phiAbsolute[K]-flowAbsolutePrevIter[K]);
                correctAbsoluteFluxBoundary(phiAbsolute[K],U[K]);
                phi[K] == phiAbsolute[K]-meshPhiHB[K];
            }

        forAll(subTimeMeshes,K) momSource[K]=Zero;
        HB.addingSource<vector>(momSource,vectorVarsGlob,0);
        hbFrame.correctSource(momSource,U,subTimeMeshes[0]);

        scalar totalVolume=0,continuitySum2=0,momentumSum2=0;
        scalar continuityMaxPerVolume=0,momentumMaxPerVolume=0;
        scalar physicalContinuityMax=0,physicalMomentumMax=0;
        bool finite=true;
        forAll(subTimeMeshes,K)
        {
            const fvMesh& mesh=subTimeMeshes[K];
            const HBPhaseTime phaseClock
            (
                runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,
                subTimeMeshes[K],meshPhiHB[K]
            );
            // Native conservative convection/stress, including patch and
            // interface sources. Never insert this temporary matrix into M.
            fvVectorMatrix stageMomentumEqn
            (
                fvm::div(phi[K],U[K])+turbulence[K].divDevReff(U[K])
            );
            const volVectorField stageH(stageMomentumEqn.H());
            const volScalarField stageA(stageMomentumEqn.A());
            const volVectorField stageGradP(fvc::grad(p[K]));
            const volScalarField stageDivQ(fvc::div(phiAbsolute[K]));
            scalarField& negativeHc=stageNegativeHc[K];
            vectorField& negativeHm=stageNegativeHm[K];
            negativeHc=-mesh.V()*stageDivQ.primitiveField();
            negativeHm=mesh.V()
                *(stageH.primitiveField()-stageA.primitiveField()*U[K].primitiveField()
                  -stageGradP.primitiveField())+momSource[K];
            forAll(negativeHc,celli)
            {
                physicalContinuityMax=max
                    (physicalContinuityMax,mag(negativeHc[celli])/mesh.V()[celli]);
                physicalMomentumMax=max
                    (physicalMomentumMax,mag(negativeHm[celli])/mesh.V()[celli]);
            }
            // -H_m = -F_m - T*deltaU. M already contains +T; subtracting it
            // here is essential when solving for a NEW TOTAL increment.
            if (useStageMomentumPseudoTime)
                negativeHm -= mesh.V()*rPseudoDeltaT[K].primitiveField()
                    *vectorSolIncrement[K].primitiveField();
            // The gauge cell row is the gauge alone: its continuity defect is
            // the closed-domain global imbalance and must not enter the
            // advective row operation either.
            // Keep Q(x) and its full Jacobian response. The integrated
            // target is fixed at (1-eta)*div(q_outer) for every inner evaluation.
            if (eta<1)
                negativeHc += stageContinuityOffset[K];
            if (gaugeReplacesContinuity && pressureNeedsReference[K]
             && pRefCells[K]>=0)
                negativeHc[pRefCells[K]]=0;
            if (useStageAdvectiveMomentum)
            {
                // Nonlinear row operation Hm <- Hm-U*Fc, before the pressure
                // gauge is added to Fc. It removes the large U*div(deltaq)
                // contribution from the momentum/pressure Jacobian. Its roots
                // are identical when Fc=0; outer acceptance ALWAYS measures
                // the original conservative momentum and face closure.
                negativeHm -= U[K].primitiveField()*negativeHc;
            }
            if (pressureNeedsReference[K] && pRefCells[K]>=0)
            {
                const label refCell=pRefCells[K];
                negativeHc[refCell] += pressureGaugeWeight[K]
                    *(pRefValues[K]-scalarVarsPrevIter[K][refCell]
                      -scalarSolIncrement[K][refCell]);
            }
            forAll(negativeHc,celli)
            {
                const scalar volume=mesh.V()[celli];
                const scalar c=mag(negativeHc[celli])/volume;
                const scalar m=mag(negativeHm[celli])/volume;
                finite=finite && std::isfinite(c) && std::isfinite(m)
                    && std::isfinite(p[K][celli]);
                for (direction cmpt=0;cmpt<vector::nComponents;++cmpt)
                    finite=finite && std::isfinite(U[K][celli][cmpt]);
                totalVolume+=volume;
                continuitySum2+=volume*sqr(c/stageContinuityScale);
                momentumSum2+=volume*sqr(m/stageMomentumScale);
                continuityMaxPerVolume=max(continuityMaxPerVolume,c);
                momentumMaxPerVolume=max(momentumMaxPerVolume,m);
            }
        }
        reduce(finite,andOp<bool>());
        reduce(totalVolume,sumOp<scalar>());
        reduce(continuitySum2,sumOp<scalar>());
        reduce(momentumSum2,sumOp<scalar>());
        reduce(continuityMaxPerVolume,maxOp<scalar>());
        reduce(momentumMaxPerVolume,maxOp<scalar>());
        reduce(physicalContinuityMax,maxOp<scalar>());
        reduce(physicalMomentumMax,maxOp<scalar>());
        stageEvaluationFinite=finite && std::isfinite(continuitySum2)
            && std::isfinite(momentumSum2) && std::isfinite(totalVolume)
            && totalVolume>0;
        if (!stageEvaluationFinite)
        {
            if (report) Info<< "IHBFullStage corrector " << corrector
                << " finite 0 scaledRms " << GREAT << endl;
            return scalar(GREAT);
        }
        const scalar cRms=Foam::sqrt(continuitySum2/totalVolume);
        const scalar mRms=Foam::sqrt(momentumSum2/totalVolume);
        const scalar stageMerit=std::hypot(cRms,mRms);
        stageEvaluationFinite=std::isfinite(cRms) && std::isfinite(mRms)
            && std::isfinite(stageMerit);
        if (!stageEvaluationFinite)
        {
            if (report) Info<< "IHBFullStage corrector " << corrector
                << " finite 0 scaledRms " << GREAT << endl;
            return scalar(GREAT);
        }
        if (report) Info<< "IHBFullStage corrector " << corrector
            << " finite 1 scaledRms " << stageMerit
            << " continuityScaledRms " << cRms
            << " momentumScaledRms " << mRms
            << " scaledMax " << max(continuityMaxPerVolume/stageContinuityScale,
                                    momentumMaxPerVolume/stageMomentumScale)
            << " continuityMaxPerVolume " << continuityMaxPerVolume
            << " momentumMaxPerVolume " << momentumMaxPerVolume
            << " physicalContinuityMax " << physicalContinuityMax
            << " physicalMomentumMax " << physicalMomentumMax << endl;
        return stageMerit;
    };

    scalar stageCurrentMerit=evaluateFullStage(0);
    if (!stageEvaluationFinite)
        FatalErrorInFunction << "Nonfinite initial complete nonlinear stage residual"
            << exit(FatalError);
    if (stageNewtonSteps>0)
    {
        #include "stageNewtonKrylov.H"
    }
    else
    {
    for (label corrector=1;corrector<=fullStageCorrectors;++corrector)
    {
        if (stageCurrentMerit<=VSMALL) break;
        forAll(subTimeMeshes,K)
        {
            stageMatrixPressure[K]=Zero;
            stageMatrixMomentum[K]=Zero;
        }
        // M includes the gauge, HB, pseudo-time and all AMI/processor actions.
        // matrixMul uses coefficients only, not the stored equation sources.
        eqSystem.matrixMul(scalarSolIncrement,vectorSolIncrement,
                           stageMatrixPressure,stageMatrixMomentum);
        PtrList<scalarField> oldPressureIncrement(nO);
        PtrList<vectorField> oldVelocityIncrement(nO);
        forAll(subTimeMeshes,K)
        {
            oldPressureIncrement.set
                (K,new scalarField(scalarSolIncrement[K].primitiveField()));
            oldVelocityIncrement.set
                (K,new vectorField(vectorSolIncrement[K].primitiveField()));
            // solveForIncr sums ALL block sources in the row. Clear them first,
            // then put the complete M*x-H in exactly one block per equation.
            forAll(subTimeMeshes,J)
            {
                eqSystem.dSByS(K,J).source()=Zero;
                eqSystem.dSByV(K,J).source()=Zero;
                eqSystem.dVByS(K,J).source()=Zero;
                eqSystem.dVByV(K,J).source()=Zero;
            }
            eqSystem.dSByS(K,K).source()=stageMatrixPressure[K]+stageNegativeHc[K];
            eqSystem.dVByV(K,K).source()=stageMatrixMomentum[K]+stageNegativeHm[K];
        }
        autoPtr<residualsIO> previousStageLinearResidual;
        if (initRes.valid())
            previousStageLinearResidual.reset(new residualsIO(initRes()));
        initRes.clear();
        initRes.set(new residualsIO(eqSystem.solveForIncr
        (
            scalarVarsPrevIter,vectorVarsPrevIter,
            scalarSolIncrement,vectorSolIncrement
        )));
        // Backtrack ONLY the new Picard correction s=xPicard-xOld. Each trial
        // remains a total increment from the same outer baseline.
        PtrList<scalarField> picardPressureIncrement(nO);
        PtrList<vectorField> picardVelocityIncrement(nO);
        forAll(subTimeMeshes,K)
        {
            picardPressureIncrement.set
                (K,new scalarField(scalarSolIncrement[K].primitiveField()));
            picardVelocityIncrement.set
                (K,new vectorField(vectorSolIncrement[K].primitiveField()));
        }
        const scalar oldStageMerit=stageCurrentMerit;
        bool stageAccepted=false;
        scalar beta=fullStageDamping;
        for (label trial=0;trial<=fullStageMaxBacktracks;++trial,beta*=0.5)
        {
            {
                stageFloatingPointGuard stageInterpolationFpeGuard;
                forAll(subTimeMeshes,K)
                {
                    // beta=1 reproduces the original successful Picard path
                    // exactly, without an unnecessary subtract/add round trip.
                    if (beta==1)
                    {
                        scalarSolIncrement[K].primitiveFieldRef()=picardPressureIncrement[K];
                        vectorSolIncrement[K].primitiveFieldRef()=picardVelocityIncrement[K];
                    }
                    else
                    {
                        scalarSolIncrement[K].primitiveFieldRef()=
                            oldPressureIncrement[K]
                            +beta*(picardPressureIncrement[K]-oldPressureIncrement[K]);
                        vectorSolIncrement[K].primitiveFieldRef()=
                            oldVelocityIncrement[K]
                            +beta*(picardVelocityIncrement[K]-oldVelocityIncrement[K]);
                    }
                }
            }
            const scalar trialMerit=evaluateFullStage(corrector);
            const scalar target=(1.0-1e-4*beta)*oldStageMerit;
            Info<< "IHBFullStageTrial corrector " << corrector
                << " trial " << trial << " beta " << beta
                << " finite " << stageEvaluationFinite
                << " scaledRms " << trialMerit << " target " << target
                << " previousScaledRms " << oldStageMerit << endl;
            if (stageEvaluationFinite && trialMerit<=target)
            {
                stageAccepted=true;
                stageCurrentMerit=trialMerit;
                Info<< "IHBFullStageAccepted corrector " << corrector
                    << " beta " << beta << " scaledRms " << stageCurrentMerit
                    << " previousScaledRms " << oldStageMerit << endl;
                break;
            }
        }
        if (!stageAccepted)
        {
            // Earlier accepted inner steps survive. Re-evaluate this saved
            // total increment to restore -H; never keep a rejected candidate.
            forAll(subTimeMeshes,K)
            {
                scalarSolIncrement[K].primitiveFieldRef()=oldPressureIncrement[K];
                vectorSolIncrement[K].primitiveFieldRef()=oldVelocityIncrement[K];
            }
            stageCurrentMerit=evaluateFullStage(corrector);
            if (!stageEvaluationFinite)
                FatalErrorInFunction << "Failed to restore the finite inner-stage baseline"
                    << exit(FatalError);
            initRes.clear();
            if (previousStageLinearResidual.valid())
                initRes.reset(new residualsIO(previousStageLinearResidual()));
            Info<< "IHBFullStageRejected corrector " << corrector
                << " restoredScaledRms " << stageCurrentMerit
                << " previousScaledRms " << oldStageMerit
                << " retainedCorrectors " << corrector-1 << endl;
            break;
        }
    }

    } // Picard alternative to stage Newton

    // No candidate becomes the physical state before the existing complete
    // (p,U,q) line search. Restore all sources and the exact saved baseline;
    // only the final total increments intentionally survive this include.
    forAll(subTimeMeshes,K)
    {
        p[K] == scalarVarsPrevIter[K];
        U[K] == vectorVarsPrevIter[K];
        phiAbsolute[K] == flowAbsolutePrevIter[K];
        phi[K] == flowRelativePrevIter[K];
        contSource[K]=stageSavedContinuitySource[K];
        momSource[K]=stageSavedMomentumSource[K];
        if (stageNewtonSteps==0)
        {
            forAll(subTimeMeshes,J)
            {
                const label block=K*nO+J;
                eqSystem.dSByS(K,J).source()=stageSavedSS[block];
                eqSystem.dSByV(K,J).source()=stageSavedSV[block];
                eqSystem.dVByS(K,J).source()=stageSavedVS[block];
                eqSystem.dVByV(K,J).source()=stageSavedVV[block];
            }
        }
    }
    // Restore flux-dependent boundary coefficients on the baseline branch too.
    forAll(subTimeMeshes,K)
    {
        const HBPhaseTime phaseClock
        (
            runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,
            subTimeMeshes[K],meshPhiHB[K]
        );
        p[K].correctBoundaryConditions();
        U[K].correctBoundaryConditions();
        MRF[K].correctBoundaryVelocity(U[K]);
    }
}
