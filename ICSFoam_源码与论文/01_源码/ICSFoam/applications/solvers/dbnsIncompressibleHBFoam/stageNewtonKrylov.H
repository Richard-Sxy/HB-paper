// Right-preconditioned, matrix-free Newton-FGMRES for the current complete H.
// Included inside fullStageDefectCorrection.H, after evaluateFullStage exists.
// x remains a TOTAL outer increment; the Newton unknown here is its correction.
{
    stageFloatingPointGuard stageNewtonArithmeticGuard;
    labelList nkOffsets(nO+1,0);
    scalar nkVolume=0;
    forAll(subTimeMeshes,K)
    {
        nkOffsets[K+1]=nkOffsets[K]+4*subTimeMeshes[K].nCells();
        forAll(subTimeMeshes[K].V(),celli) nkVolume+=subTimeMeshes[K].V()[celli];
    }
    reduce(nkVolume,sumOp<scalar>());
    if (!std::isfinite(nkVolume) || nkVolume<=0)
        FatalErrorInFunction << "Invalid Newton/Krylov volume" << exit(FatalError);
    const label nkSize=nkOffsets[nO];
    const scalar nkPressureScale=sqr(flowVelocityScale);
    const scalar nkVelocityScale=flowVelocityScale;
    scalarField nkWeights(nkSize,Zero);
    forAll(subTimeMeshes,K)
        forAll(subTimeMeshes[K].V(),celli)
            for (label cmpt=0;cmpt<4;++cmpt)
                nkWeights[nkOffsets[K]+4*celli+cmpt]=subTimeMeshes[K].V()[celli]/nkVolume;

    const auto nkDot = [&](const scalarField& a,const scalarField& b)
    {
        scalar value=0;
        forAll(a,i) value+=nkWeights[i]*a[i]*b[i];
        reduce(value,sumOp<scalar>());
        return value;
    };
    const auto nkNorm = [&](const scalarField& a)
    {
        const scalar square=nkDot(a,a);
        return (!std::isfinite(square) || square<0)
            ? scalar(-1) : Foam::sqrt(square);
    };
    const auto nkFinite = [&](const scalarField& a)
    {
        bool finite=true;
        forAll(a,i) finite=finite && std::isfinite(a[i]);
        reduce(finite,andOp<bool>());
        return finite;
    };
    const auto nkPackNegativeH = [&]()
    {
        scalarField result(nkSize,Zero);
        forAll(subTimeMeshes,K)
            forAll(stageNegativeHc[K],celli)
            {
                const label start=nkOffsets[K]+4*celli;
                const scalar volume=subTimeMeshes[K].V()[celli];
                result[start]=stageNegativeHc[K][celli]/(volume*stageContinuityScale);
                for (direction cmpt=0;cmpt<3;++cmpt)
                    result[start+1+cmpt]=stageNegativeHm[K][celli].component(cmpt)
                        /(volume*stageMomentumScale);
            }
        return result;
    };

    // Only the public solver/preconditioner API is used. Its assembled M and
    // cached local inverse stay fixed for this entire nonlinear stage.
    const dictionary& nkFlowDict=eqSystem.mesh().solutionDict().subDict("flowSolver");
    const word nkSolverName(nkFlowDict.get<word>("solver"));
    dictionary nkPreconditionerControls(nkFlowDict.subDict(nkSolverName));
    nkPreconditionerControls.set("preconditioner",word("ihbBlockJacobi"));
    autoPtr<coupledMatrix::solver> nkOwner;
    autoPtr<coupledMatrix::preconditioner> nkPreconditioner;
    if (stagePreconditioner=="blockJacobi")
    {
        nkOwner=coupledMatrix::solver::New(nkFlowDict,eqSystem);
        nkPreconditioner=coupledMatrix::preconditioner::New
            (nkOwner(),nkPreconditionerControls);
    }
    PtrList<scalarField> nkScalarWork(nO);
    PtrList<vectorField> nkVectorWork(nO);
    forAll(subTimeMeshes,K)
    {
        nkScalarWork.set(K,new scalarField(subTimeMeshes[K].nCells(),Zero));
        nkVectorWork.set(K,new vectorField(subTimeMeshes[K].nCells(),Zero));
    }
    #include "stagePressureSchur.H"

    const auto nkApplyPreconditioner = [&](const scalarField& input,scalarField& output)
    {
        // P_hat^-1 = X^-1 * P^-1 * R^-1, where R scales integrated residuals
        // by division through V*Cref or V*Mref, while X uses Uref^2 or Uref.
        forAll(subTimeMeshes,K)
            forAll(nkScalarWork[K],celli)
            {
                const label start=nkOffsets[K]+4*celli;
                const scalar volume=subTimeMeshes[K].V()[celli];
                nkScalarWork[K][celli]=volume*stageContinuityScale*input[start];
                for (direction cmpt=0;cmpt<3;++cmpt)
                    nkVectorWork[K][celli].component(cmpt)=
                        volume*stageMomentumScale*input[start+1+cmpt];
            }
        if (stagePreconditioner=="pressureSchur")
        {
            if (!nkApplyPressureSchur()) return false;
        }
        else nkPreconditioner->precondition(nkScalarWork,nkVectorWork);
        output.setSize(nkSize);
        forAll(subTimeMeshes,K)
            forAll(nkScalarWork[K],celli)
            {
                const label start=nkOffsets[K]+4*celli;
                output[start]=nkScalarWork[K][celli]/nkPressureScale;
                for (direction cmpt=0;cmpt<3;++cmpt)
                    output[start+1+cmpt]=nkVectorWork[K][celli].component(cmpt)/nkVelocityScale;
            }
        return nkFinite(output);
    };

    const scalar nkInitialMerit=stageCurrentMerit;
    const scalar nkNonlinearTarget=stageNewtonRelTol*nkInitialMerit;
    Info<< "IHBStageNewtonStart scaledRms " << nkInitialMerit
        << " target " << nkNonlinearTarget
        << " directions " << stageKrylovDirections << " cycles " << stageKrylovCycles
        << " preconditioner " << stagePreconditioner
        << " side right method FGMRES" << endl;

    for (label newton=1;newton<=stageNewtonSteps;++newton)
    {
        if (stageCurrentMerit<=nkNonlinearTarget || stageCurrentMerit<=VSMALL) break;
        PtrList<scalarField> nkBasePressure(nO),nkBaseNegativeHc(nO);
        PtrList<vectorField> nkBaseVelocity(nO),nkBaseNegativeHm(nO);
        scalar stateMaximum=0;
        forAll(subTimeMeshes,K)
        {
            nkBasePressure.set(K,new scalarField(scalarSolIncrement[K].primitiveField()));
            nkBaseVelocity.set(K,new vectorField(vectorSolIncrement[K].primitiveField()));
            nkBaseNegativeHc.set(K,new scalarField(stageNegativeHc[K]));
            nkBaseNegativeHm.set(K,new vectorField(stageNegativeHm[K]));
            forAll(nkBasePressure[K],celli)
            {
                stateMaximum=max(stateMaximum,mag
                    ((scalarVarsPrevIter[K][celli]+nkBasePressure[K][celli])/nkPressureScale));
                for (direction cmpt=0;cmpt<3;++cmpt)
                    stateMaximum=max(stateMaximum,mag
                    (
                        (vectorVarsPrevIter[K][celli].component(cmpt)
                        +nkBaseVelocity[K][celli].component(cmpt))/nkVelocityScale
                    ));
            }
        }
        reduce(stateMaximum,maxOp<scalar>());
        const scalar nkOldMerit=stageCurrentMerit;
        const scalarField rhs(nkPackNegativeH()); // R*(-H) at this Newton base.
        const scalar rhsNorm=nkNorm(rhs);
        if (!nkFinite(rhs) || !std::isfinite(stateMaximum) || rhsNorm<0)
            FatalErrorInFunction << "Nonfinite stage Newton baseline" << exit(FatalError);
        const scalar packedNormRelativeError=mag(rhsNorm-nkOldMerit)
            /max(max(rhsNorm,nkOldMerit),scalar(VSMALL));
        if (!std::isfinite(packedNormRelativeError) || packedNormRelativeError>1e-10)
            FatalErrorInFunction
                << "Packed stage residual norm differs from complete H merit: "
                << rhsNorm << " versus " << nkOldMerit
                << ", relative error " << packedNormRelativeError << exit(FatalError);
        if (newton==1)
            Info<< "IHBStageNewtonNormCheck packedScaledRms " << rhsNorm
                << " completeScaledRms " << nkOldMerit
                << " relativeError " << packedNormRelativeError << endl;
        if (rhsNorm<=VSMALL) break;
        const scalar linearTarget=stageKrylovRelTol*rhsNorm;

        const auto nkRestoreBase = [&]()
        {
            forAll(subTimeMeshes,K)
            {
                scalarSolIncrement[K].primitiveFieldRef()=nkBasePressure[K];
                vectorSolIncrement[K].primitiveFieldRef()=nkBaseVelocity[K];
                stageNegativeHc[K]=nkBaseNegativeHc[K];
                stageNegativeHm[K]=nkBaseNegativeHm[K];
            }
            stageEvaluationFinite=true;
        };
        const auto nkSetProbe = [&](const scalarField& searchDirection,const scalar factor)
        {
            forAll(subTimeMeshes,K)
            {
                scalarField& pressureIncrement=scalarSolIncrement[K].primitiveFieldRef();
                vectorField& velocityIncrement=vectorSolIncrement[K].primitiveFieldRef();
                forAll(nkBasePressure[K],celli)
                {
                    const label start=nkOffsets[K]+4*celli;
                    pressureIncrement[celli]=nkBasePressure[K][celli]
                        +factor*nkPressureScale*searchDirection[start];
                    for (direction cmpt=0;cmpt<3;++cmpt)
                        velocityIncrement[celli].component(cmpt)=
                            nkBaseVelocity[K][celli].component(cmpt)
                           +factor*nkVelocityScale*searchDirection[start+1+cmpt];
                }
            }
        };
        label fdEvaluations=0;
        const auto nkApplyJacobian = [&](const scalarField& searchDirection,scalarField& result)
        {
            result.setSize(nkSize);
            result=Zero;
            if (!nkFinite(searchDirection)) return false;
            scalar maximum=0;
            forAll(searchDirection,i) maximum=max(maximum,mag(searchDirection[i]));
            reduce(maximum,maxOp<scalar>());
            if (maximum<=VSMALL) return true;
            const scalar h=stageFiniteDifferenceEpsilon*(1+stateMaximum)/maximum;
            if (!std::isfinite(h) || h<=0) return false;
            nkSetProbe(searchDirection,h);
            evaluateFullStage(newton,false);
            ++fdEvaluations;
            bool valid=stageEvaluationFinite;
            if (valid)
            {
                const scalarField perturbed(nkPackNegativeH());
                // stageNegativeH is -H: J*v = (-H_base - (-H_probe))/h.
                result=(rhs-perturbed)/h;
                valid=nkFinite(result);
            }
            nkRestoreBase();
            return valid;
        };

        scalarField totalStep(nkSize,Zero),bestStep(nkSize,Zero),residual(rhs);
        scalar bestLinearNorm=rhsNorm;
        scalar bestDirectionalDescent=0;
        bool krylovFinite=true,linearConverged=false;
        for (label cycle=0;cycle<stageKrylovCycles;++cycle)
        {
            const scalar beta=nkNorm(residual);
            if (beta<0) { krylovFinite=false; break; }
            if (beta<=linearTarget) { linearConverged=true; break; }
            const label m=stageKrylovDirections;
            List<scalarField> basis(m+1),preconditioned(m);
            List<scalarField> hessenberg(m+1);
            for (label row=0;row<=m;++row)
                hessenberg[row]=scalarField(m,Zero);
            scalarField cosine(m,Zero),sine(m,Zero),g(m+1,Zero);
            basis[0]=residual/beta;
            g[0]=beta;
            label used=0;
            bool happyBreakdown=false;
            scalar estimatedNorm=beta;
            for (label column=0;column<m;++column)
            {
                if (!nkApplyPreconditioner(basis[column],preconditioned[column]))
                { krylovFinite=false; break; }
                scalarField work(nkSize,Zero);
                if (!nkApplyJacobian(preconditioned[column],work))
                { krylovFinite=false; break; }
                const scalar originalNorm=nkNorm(work);
                if (originalNorm<0) { krylovFinite=false; break; }
                // Twice-reorthogonalised modified Gram-Schmidt in the SAME
                // global volume-weighted metric as the nonlinear stage merit.
                for (label sweep=0;sweep<2 && krylovFinite;++sweep)
                    for (label row=0;row<=column;++row)
                    {
                        const scalar projection=nkDot(work,basis[row]);
                        if (!std::isfinite(projection))
                        { krylovFinite=false; break; }
                        hessenberg[row][column]+=projection;
                        work-=projection*basis[row];
                    }
                if (!krylovFinite) break;
                const scalar nextNorm=nkNorm(work);
                if (nextNorm<0) { krylovFinite=false; break; }
                hessenberg[column+1][column]=nextNorm;
                happyBreakdown=nextNorm<=1e-12*max(originalNorm,scalar(VSMALL));
                for (label row=0;row<column;++row)
                {
                    const scalar top=hessenberg[row][column];
                    const scalar bottom=hessenberg[row+1][column];
                    hessenberg[row][column]=cosine[row]*top+sine[row]*bottom;
                    hessenberg[row+1][column]=-sine[row]*top+cosine[row]*bottom;
                }
                const scalar diagonal=hessenberg[column][column];
                const scalar subdiagonal=hessenberg[column+1][column];
                const scalar rotationNorm=std::hypot(diagonal,subdiagonal);
                if (!std::isfinite(rotationNorm)) { krylovFinite=false; break; }
                if (rotationNorm<=VSMALL) break; // Singular direction, no division.
                cosine[column]=diagonal/rotationNorm;
                sine[column]=subdiagonal/rotationNorm;
                hessenberg[column][column]=rotationNorm;
                hessenberg[column+1][column]=0;
                const scalar previousG=g[column];
                g[column]=cosine[column]*previousG;
                g[column+1]=-sine[column]*previousG;
                estimatedNorm=mag(g[column+1]);
                used=column+1;
                if (happyBreakdown || estimatedNorm<=linearTarget) break;
                basis[column+1]=work/nextNorm;
            }
            if (!krylovFinite || !used) break;
            scalarField coefficients(used,Zero);
            for (label row=used-1;row>=0;--row)
            {
                scalar value=g[row];
                for (label col=row+1;col<used;++col)
                    value-=hessenberg[row][col]*coefficients[col];
                const scalar diagonal=hessenberg[row][row];
                if (!std::isfinite(value) || !std::isfinite(diagonal) || mag(diagonal)<=VSMALL)
                { krylovFinite=false; break; }
                coefficients[row]=value/diagonal;
            }
            if (!krylovFinite) break;
            for (label column=0;column<used;++column)
                totalStep+=coefficients[column]*preconditioned[column];
            if (!nkFinite(totalStep)) { krylovFinite=false; break; }
            scalarField trueAction(nkSize,Zero);
            if (!nkApplyJacobian(totalStep,trueAction)) { krylovFinite=false; break; }
            residual=rhs-trueAction;
            const scalar trueNorm=nkNorm(residual);
            if (trueNorm<0) { krylovFinite=false; break; }
            if (trueNorm<bestLinearNorm)
            {
                // d||H||/d beta = -<rhs,J*step>/||rhs||. Use the independently
                // recomputed FD action, not the Hessenberg residual estimate.
                const scalar descent=nkDot(rhs,trueAction)/rhsNorm;
                if (!std::isfinite(descent) || descent<=0)
                { krylovFinite=false; break; }
                bestLinearNorm=trueNorm;
                bestDirectionalDescent=descent;
                bestStep=totalStep;
            }
            Info<< "IHBStageKrylov newton " << newton << " cycle " << cycle+1
                << " directions " << used << " happyBreakdown " << happyBreakdown
                << " estimatedRelative " << estimatedNorm/rhsNorm
                << " trueFdRelative " << trueNorm/rhsNorm
                << " targetRelative " << stageKrylovRelTol << endl;
            if (trueNorm<=linearTarget) { linearConverged=true; break; }
            // Restart from the best verified linear correction, not a worse
            // finite-difference/noise-affected last cycle.
            if (trueNorm>bestLinearNorm)
            {
                totalStep=bestStep;
                if (!nkApplyJacobian(totalStep,trueAction)) { krylovFinite=false; break; }
                residual=rhs-trueAction;
            }
        }

        // Restore actual fields/patch state as well as the cached increment/H.
        // Every probe itself was rebuilt from the same outer baseline.
        nkRestoreBase();
        const scalar restoredMerit=evaluateFullStage(newton,false);
        if (!stageEvaluationFinite
         || mag(restoredMerit-nkOldMerit)>1e-9*max(scalar(1),nkOldMerit))
            FatalErrorInFunction << "Newton FD probes did not reproduce their baseline"
                << exit(FatalError);
        stageCurrentMerit=restoredMerit;
        if (!krylovFinite || bestLinearNorm>=rhsNorm)
        {
            Info<< "IHBStageNewtonRejected newton " << newton
                << " reason invalidOrNonreducingLinearStep"
                << " finite " << krylovFinite << " fdEvaluations " << fdEvaluations
                << " retainedScaledRms " << stageCurrentMerit << endl;
            break;
        }

        bool accepted=false;
        scalar beta=fullStageDamping;
        for (label trial=0;trial<=fullStageMaxBacktracks;++trial,beta*=0.5)
        {
            nkSetProbe(bestStep,beta);
            const scalar trialMerit=evaluateFullStage(newton);
            const scalar target=nkOldMerit
                -stageNewtonArmijo*beta*bestDirectionalDescent;
            Info<< "IHBStageNewtonTrial newton " << newton << " trial " << trial
                << " beta " << beta << " finite " << stageEvaluationFinite
                << " scaledRms " << trialMerit << " target " << target
                << " directionalDescent " << bestDirectionalDescent
                << " linearConverged " << linearConverged
                << " trueFdRelative " << bestLinearNorm/rhsNorm << endl;
            if (stageEvaluationFinite && trialMerit<=target)
            {
                accepted=true;
                stageCurrentMerit=trialMerit;
                Info<< "IHBStageNewtonAccepted newton " << newton
                    << " beta " << beta << " scaledRms " << stageCurrentMerit
                    << " previousScaledRms " << nkOldMerit
                    << " fdEvaluations " << fdEvaluations << endl;
                break;
            }
        }
        if (!accepted)
        {
            nkRestoreBase();
            stageCurrentMerit=evaluateFullStage(newton);
            if (!stageEvaluationFinite)
                FatalErrorInFunction << "Failed to restore rejected Newton stage"
                    << exit(FatalError);
            Info<< "IHBStageNewtonRejected newton " << newton
                << " reason nonlinearBacktracking"
                << " retainedScaledRms " << stageCurrentMerit << endl;
            break;
        }
    }
}
