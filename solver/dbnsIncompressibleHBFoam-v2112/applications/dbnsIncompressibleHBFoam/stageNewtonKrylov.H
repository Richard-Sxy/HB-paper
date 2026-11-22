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
    const auto nkMomentumDot = [&](const scalarField& a,const scalarField& b)
    {
        scalar value=0;
        for (label i=0;i<nkSize;i+=4)
            for (label cmpt=1;cmpt<4;++cmpt)
                value+=nkWeights[i+cmpt]*a[i+cmpt]*b[i+cmpt];
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
    const bool matchPhysicalDiagonal=stageSchurMatchMomentumDiagonal
        && !useStageMomentumPseudoTime;
    PtrList<tensorField> savedVelocityDiagonal(nO);
    if (matchPhysicalDiagonal)
    {
        // Use the diagonal of the physical velocity block for this
        // preconditioner. The exact stage residual and Jv stay unchanged.
        forAll(subTimeMeshes,K)
        {
            tensorField& block=eqSystem.dVByV(K,K).diag();
            savedVelocityDiagonal.set(K,new tensorField(block));
            block-=subTimeMeshes[K].V()*rPseudoDeltaT[K].primitiveField()*tensor::I;
        }
        Info<< "IHBStageSchurMomentumDiagonal physical" << endl;
    }
    #include "stagePressureSchur.H"

    enum class VelocityProbe { production, assembled, frozenFlux, coupled };
    VelocityProbe velocityProbe=VelocityProbe::production;
    const auto velocityOnly = [&]()
    {
        return velocityProbe==VelocityProbe::assembled || velocityProbe==VelocityProbe::frozenFlux;
    };

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
        if (velocityOnly())
        {
            if (!nkSchurApplyVelocityInverse(nkVectorWork)) return false;
        }
        else if (stagePreconditioner=="pressureSchur")
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

    scalarField acceptedMomentumDirection;
    for (label newton=1;;++newton)
    {
        const bool continuityCorrection=stageCoupledContinuityCorrection
            && (newton>stageNewtonSteps || stageCurrentMerit<=nkNonlinearTarget);
        if (stageCurrentMerit<=VSMALL
         || (!continuityCorrection
          && (newton>stageNewtonSteps || stageCurrentMerit<=nkNonlinearTarget)))
            break;
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
        const bool audit=runTime.timeIndex()==2 && nO>2
            && physicalControls.lookupOrDefault<bool>("stageFoilAudit",false);
        const label hotPhase=2;
        const label hotRank=3;
        const label hotCell=4;
        const bool localHotspot=audit && Pstream::myProcNo()==hotRank
            && subTimeMeshes[hotPhase].nCells()>hotCell;
        const auto hotspotGradient = [&]()
        {
            const HBPhaseTime clock
            (
                runTime,snapshots[hotPhase]+hbMotionDeltaT,hbMotionDeltaT,
                subTimeMeshes[hotPhase],meshPhiHB[hotPhase]
            );
            const volVectorField gradient(fvc::grad(p[hotPhase]));
            return localHotspot ? gradient[hotCell] : vector::zero;
        };
        const auto auditFlux = [&](const char* state, const surfaceScalarField& q,
                                   const label K)
        {
            const fvMesh& mesh=subTimeMeshes[K];
            const label out=mesh.boundaryMesh().findPatchID("AMI_OUT");
            const label in=mesh.boundaryMesh().findPatchID("AMI_IN");
            if (min(out,in)<0)
                FatalErrorInFunction << "Missing AMI_OUT/AMI_IN in foil audit"
                    << exit(FatalError);
            const scalar outNet=gSum(q.boundaryField()[out]);
            const scalar inNet=gSum(q.boundaryField()[in]);
            const volScalarField divQ(fvc::div(q));
            const scalar domainNet=gSum
            (
                scalarField(mesh.V()*divQ.primitiveField())
            );
            const label ref=pressureNeedsReference[K] ? pRefCells[K] : -1;
            if (ref>=0)
                Pout<< "IHBRefContinuity state " << state << " phase " << K
                    << " rank " << Pstream::myProcNo() << " cell " << ref
                    << " divQ " << divQ[ref]
                    << " integrated " << mesh.V()[ref]*divQ[ref] << endl;
            Info<< "IHBAMINet state " << state << " phase " << K
                << " AMI_OUT " << outNet << " AMI_IN " << inNet
                << " amiSum " << outNet+inNet << " domainNet " << domainNet
                << " nonAmiNet " << domainNet-(outNet+inNet)
                << " pressureReference " << pressureNeedsReference[K] << endl;
        };
        if (audit && newton==1)
        {
            forAll(subTimeMeshes,K)
            {
                auditFlux("stored",flowAbsolutePrevIter[K],K);
                auditFlux("gmres",phiAbsolute[K],K);
                scalarSolIncrement[K].primitiveFieldRef()=Zero;
                vectorSolIncrement[K].primitiveFieldRef()=Zero;
            }
            evaluateFullStage(newton,false);
            if (!stageEvaluationFinite)
                FatalErrorInFunction << "Invalid zero-increment foil audit state"
                    << exit(FatalError);
            forAll(subTimeMeshes,K)
                auditFlux("outerClosure",phiAbsolute[K],K);
            nkRestoreBase();
            evaluateFullStage(newton,false);
            if (!stageEvaluationFinite
             || nkNorm(scalarField(rhs-nkPackNegativeH()))>1e-12*rhsNorm)
                FatalErrorInFunction << "Foil flux audit did not restore Newton base"
                    << exit(FatalError);
        }
        const vector baseGradient=audit ? hotspotGradient() : vector::zero;
        label fdEvaluations=0;
        const bool auditConvection=stageAuditConvectionFlux && !continuityCorrection;
        const bool applyConvectionFeedback=stageSchurConvectionFeedback
            && !continuityCorrection;
        const bool coupledVelocityFeedback=stageSchurVelocityConvectionFeedback
            && stagePreconditioner=="pressureSchur" && !continuityCorrection
            && !useStageAdvectiveMomentum;
        PtrList<surfaceScalarField> newtonBaseFlux(nO);
        if (auditConvection || applyConvectionFeedback || coupledVelocityFeedback)
            forAll(phi,K)
                newtonBaseFlux.set
                (
                    K,new surfaceScalarField
                    (
                        IOobject("newtonBaseRelativeFlux",runTime.timeName(),
                            subTimeMeshes[K],IOobject::NO_READ,IOobject::NO_WRITE,false),
                        phi[K]
                    )
                );

        const bool auditVelocity=stageAuditVelocityBlock && newton==1
            && !continuityCorrection && !useStageMomentumPseudoTime
            && !useStageAdvectiveMomentum;
        PtrList<surfaceScalarField> velocityProbeFlux(nO);
        PtrList<vectorField> assembledVelocityResponse(nO);
        PtrList<blockFvMatrix<vector,tensor>> velocityAuditMatchedBlocks(nO);
        const bool matchedAudit=auditVelocity && stageVelocityAuditMatchFlux;
        PtrList<blockFvMatrix<vector,tensor>> auditConvectionBlocks(nO),auditStressBlocks(nO);
        PtrList<vectorField> auditBaseConvection(nO),auditBaseStress(nO),auditBaseTemporal(nO);
        PtrList<tensorField> auditHBDiagonal(nO);
        List<scalarField> velocityBlockMismatch(3);
        bool captureVelocityMismatch=false;
        bool captureConvectionMismatch=false;
        bool applyVelocityReconstructionFeedback=false;
        const auto velocityTerm = [](const fvVectorMatrix& equation,const volVectorField& field)
        {
            const volVectorField defect(equation.A()*field-equation.H());
            return vectorField(defect.primitiveField());
        };
        if (auditVelocity)
            forAll(subTimeMeshes,K)
            {
                velocityProbeFlux.set
                (
                    K,new surfaceScalarField
                    (
                        IOobject("velocityAuditAbsoluteFlux",runTime.timeName(),
                            subTimeMeshes[K],IOobject::NO_READ,IOobject::NO_WRITE,false),
                        phiAbsolute[K]
                    )
                );
                fvMesh& mesh=subTimeMeshes[K];
                assembledVelocityResponse.set(K,new vectorField(mesh.nCells(),Zero));
                if (matchedAudit)
                {
                    const HBPhaseTime phaseClock
                    (
                        runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,mesh,meshPhiHB[K]
                    );
                    velocityAuditMatchedBlocks.set(K,new blockFvMatrix<vector,tensor>(mesh));
                    fvVectorMatrix equation
                    (
                        fvm::div(phi[K],U[K])+turbulence[K].divDevReff(U[K])
                    );
                    auto& block=velocityAuditMatchedBlocks[K];
                    block.insertEquation(equation);
                    blockFvMatrix<vector,tensor> hbDiagonal(mesh);
                    hbDiagonal.diag()=Zero;
                    HB.addBlock(hbDiagonal,K,K);
                    hbFrame.correctBlock(hbDiagonal.diag(),K,K,mesh);
                    block.diag()+=hbDiagonal.diag();

                    fvVectorMatrix convection(fvm::div(phi[K],U[K]));
                    fvVectorMatrix stress(turbulence[K].divDevReff(U[K]));
                    auditBaseConvection.set(K,new vectorField(velocityTerm(convection,U[K])));
                    auditBaseStress.set(K,new vectorField(velocityTerm(stress,U[K])));
                    auditBaseTemporal.set(K,new vectorField(-momSource[K]/mesh.V()));
                    auditHBDiagonal.set(K,new tensorField(hbDiagonal.diag()));
                    auditConvectionBlocks.set(K,new blockFvMatrix<vector,tensor>(mesh));
                    auditStressBlocks.set(K,new blockFvMatrix<vector,tensor>(mesh));
                    auditConvectionBlocks[K].insertEquation(convection);
                    auditStressBlocks[K].insertEquation(stress);
                }
            }
        if (coupledVelocityFeedback && !matchedAudit)
            forAll(subTimeMeshes,K)
            {
                fvMesh& mesh=subTimeMeshes[K];
                const HBPhaseTime phaseClock
                (
                    runTime,snapshots[K]+hbMotionDeltaT,hbMotionDeltaT,mesh,meshPhiHB[K]
                );
                fvVectorMatrix convection(fvm::div(phi[K],U[K]));
                auditBaseConvection.set(K,new vectorField(velocityTerm(convection,U[K])));
                auditConvectionBlocks.set(K,new blockFvMatrix<vector,tensor>(mesh));
                auditConvectionBlocks[K].insertEquation(convection);
            }
        scalar maxFrozenMassAction=0;

        const auto nkApplyJacobian = [&]
        (
            const scalarField& searchDirection,
            scalarField& result,
            scalarField* convectiveFluxAction=nullptr
        )
        {
            result.setSize(nkSize);
            result=Zero;
            if (convectiveFluxAction)
            {
                convectiveFluxAction->setSize(nkSize);
                *convectiveFluxAction=Zero;
            }
            if (!nkFinite(searchDirection)) return false;
            if (velocityProbe==VelocityProbe::assembled)
            {
                forAll(subTimeMeshes,K)
                    forAll(nkVectorWork[K],celli)
                    {
                        const label i=nkOffsets[K]+4*celli;
                        for (direction c=0;c<3;++c)
                            nkVectorWork[K][celli].component(c)=nkVelocityScale*searchDirection[i+c+1];
                    }
                nkSchurVelocityAction
                (
                    nkVectorWork,assembledVelocityResponse,
                    matchedAudit ? &velocityAuditMatchedBlocks : nullptr
                );
                forAll(subTimeMeshes,K)
                    forAll(nkVectorWork[K],celli)
                    {
                        const scalar volume=subTimeMeshes[K].V()[celli];
                        vector action=assembledVelocityResponse[K][celli]/volume;
                        if (!matchPhysicalDiagonal && !matchedAudit)
                            action-=rPseudoDeltaT[K][celli]*nkVectorWork[K][celli];
                        const label i=nkOffsets[K]+4*celli;
                        for (direction c=0;c<3;++c)
                            result[i+c+1]=action.component(c)/stageMomentumScale;
                    }
                return nkFinite(result);
            }
            scalar maximum=0;
            forAll(searchDirection,i) maximum=max(maximum,mag(searchDirection[i]));
            reduce(maximum,maxOp<scalar>());
            if (maximum<=VSMALL) return true;
            const scalar h=stageFiniteDifferenceEpsilon*(1+stateMaximum)/maximum;
            if (!std::isfinite(h) || h<=0) return false;
            nkSetProbe(searchDirection,h);
            evaluateFullStage
            (
                newton,false,
                velocityProbe==VelocityProbe::frozenFlux ? &velocityProbeFlux : nullptr
            );
            ++fdEvaluations;
            bool valid=stageEvaluationFinite;
            if (valid)
            {
                const scalarField perturbed(nkPackNegativeH());
                // stageNegativeH is -H: J*v = (-H_base - (-H_probe))/h.
                result=(rhs-perturbed)/h;
                if (velocityProbe==VelocityProbe::frozenFlux)
                    for (label i=0;i<nkSize;i+=4)
                    {
                        maxFrozenMassAction=max(maxFrozenMassAction,mag(result[i]));
                        result[i]=0;
                    }
                valid=nkFinite(result);
                if (valid && (captureVelocityMismatch || captureConvectionMismatch))
                {
                    #include "stageVelocityDefectTerms.H"
                }
                if (valid && convectiveFluxAction)
                {
                    #include "stageConvectionFluxAction.H"
                    valid=nkFinite(*convectiveFluxAction);
                }
            }
            nkRestoreBase();
            return valid;
        };


        const auto solveFrozenRhs = [&]
        (
            const scalarField& rhs, const scalar rhsNorm,
            const scalar linearTarget, const char* rhsKind,
            scalarField& bestStep, scalar& bestLinearNorm,
            scalar& bestDirectionalDescent, bool& krylovFinite,
            bool& linearConverged
        )
        {
            scalarField totalStep(nkSize,Zero),residual(rhs);
            bestLinearNorm=rhsNorm;
            bestDirectionalDescent=0;
            krylovFinite=true;
            linearConverged=false;
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
                    scalarField work(nkSize,Zero),fluxAction;
                    const bool productionProbe=velocityProbe==VelocityProbe::production;
                    const bool firstAudit=auditConvection && productionProbe
                        && cycle==0 && column==0;
                    const bool feedback=applyConvectionFeedback && productionProbe;
                    const bool velocityFeedback=
                        (applyVelocityReconstructionFeedback
                         && velocityProbe==VelocityProbe::frozenFlux)
                        || (coupledVelocityFeedback && productionProbe);
                    captureConvectionMismatch=velocityFeedback;
                    const bool actionValid=nkApplyJacobian
                    (
                        preconditioned[column],work,
                        feedback || firstAudit ? &fluxAction : nullptr
                    );
                    captureConvectionMismatch=false;
                    if (!actionValid)
                    { krylovFinite=false; break; }

                    if (firstAudit)
                    {
                        const scalarField frozenAction(work-fluxAction);
                        Info<< "IHBSchurConvectionFirst fullMomentumSlope "
                            << nkMomentumDot(basis[column],work)
                            << " frozenMomentumSlope "
                            << nkMomentumDot(basis[column],frozenAction)
                            << " fluxMomentumSlope "
                            << nkMomentumDot(basis[column],fluxAction)
                            << " fullActionNorm " << nkNorm(work)
                            << " fluxActionNorm " << nkNorm(fluxAction) << endl;
                    }

                    if (feedback || velocityFeedback)
                    {
                        // Approximate the inverse of the missing momentum
                        // block action with one more application of the same inverse.
                        const scalarField& missingAction=velocityFeedback
                            ? velocityBlockMismatch[0] : fluxAction;
                        const scalarField feedbackRhs(-missingAction);
                        scalarField correction,correctionAction(nkSize,Zero);
                        if (nkApplyPreconditioner(feedbackRhs,correction)
                         && nkApplyJacobian(correction,correctionAction))
                        {
                            const scalarField defect(basis[column]-work);
                            const scalar before=nkNorm(defect);
                            const scalar actionSquare=nkDot(correctionAction,correctionAction);
                            if (before>=0 && actionSquare>VSMALL && std::isfinite(actionSquare))
                            {
                                const scalar scale=nkDot(defect,correctionAction)/actionSquare;
                                if (std::isfinite(scale) && mag(scale)>SMALL)
                                {
                                    const scalarField candidate
                                        (preconditioned[column]+scale*correction);
                                    scalarField candidateAction(nkSize,Zero);
                                    if (nkApplyJacobian(candidate,candidateAction))
                                    {
                                        const scalar after=nkNorm
                                            (scalarField(basis[column]-candidateAction));
                                        const bool useCorrection=after>=0 && after<before;
                                        if (useCorrection)
                                        {
                                            preconditioned[column]=candidate;
                                            work=candidateAction;
                                        }
                                        Info<< (velocityFeedback
                                            ? (productionProbe
                                                ? "IHBCoupledVelocityReconstructionFeedback"
                                                : "IHBVelocityReconstructionFeedback")
                                            : "IHBSchurConvectionFeedback")
                                            << " cycle " << cycle+1
                                            << " column " << column+1 << " multiplier " << scale
                                            << " defectBefore " << before << " defectAfter " << after
                                            << " used " << useCorrection << endl;
                                    }
                                }
                            }
                        }
                    }
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
                scalarField trueAction(nkSize,Zero),cycleFluxAction;
                const bool auditProduction=auditConvection
                    && velocityProbe==VelocityProbe::production;
                if (!nkApplyJacobian
                    (totalStep,trueAction,auditProduction ? &cycleFluxAction : nullptr))
                { krylovFinite=false; break; }
                residual=rhs-trueAction;
                if (auditProduction)
                {
                    const scalar momentumRhsSquare=nkMomentumDot(rhs,rhs);
                    const scalarField frozenResidual(residual+cycleFluxAction);
                    Info<< "IHBKrylovMomentumAudit cycle " << cycle+1
                        << " fullRelative "
                        << Foam::sqrt(nkMomentumDot(residual,residual)
                            /max(momentumRhsSquare,scalar(VSMALL)))
                        << " frozenFluxRelative "
                        << Foam::sqrt(nkMomentumDot(frozenResidual,frozenResidual)
                            /max(momentumRhsSquare,scalar(VSMALL)))
                        << " fluxSlope "
                        << nkMomentumDot(rhs,cycleFluxAction)
                            /max(momentumRhsSquare,scalar(VSMALL))
                        << " fluxResponseRelative "
                        << Foam::sqrt(nkMomentumDot(cycleFluxAction,cycleFluxAction)
                            /max(momentumRhsSquare,scalar(VSMALL))) << endl;
                }
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
                Info<< "IHBStageKrylov newton " << newton
                    << " rhsKind " << rhsKind << " cycle " << cycle+1
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

        };
        if (auditVelocity)
        {
            #include "stageVelocityBlockAudit.H"
        }

        scalarField continuityRhs;
        if (continuityCorrection)
        {
            continuityRhs=rhs;
            // Retain all Jacobian rows, including the full HB momentum
            // response. Drive only the continuity defect in this correction.
            for (label i=0;i<nkSize;i+=4)
                for (label component=1;component<4;++component)
                    continuityRhs[i+component]=0;
        }
        const scalarField& correctionRhs=continuityCorrection ? continuityRhs : rhs;
        const scalar correctionRhsNorm=nkNorm(correctionRhs);
        if (correctionRhsNorm<=VSMALL) break;
        if (continuityCorrection && stageAuditContinuityCorrection)
        {
            scalar scalarMaximum=0,momentumMaximum=0;
            for (label i=0;i<nkSize;i+=4)
            {
                scalarMaximum=max(scalarMaximum,mag(continuityRhs[i]));
                for (label cmpt=1;cmpt<4;++cmpt)
                    momentumMaximum=max(momentumMaximum,mag(continuityRhs[i+cmpt]));
            }
            reduce(scalarMaximum,maxOp<scalar>());
            reduce(momentumMaximum,maxOp<scalar>());
            Info<< "IHBContinuityRhsAudit rms " << correctionRhsNorm
                << " max " << scalarMaximum << " maxMomentum " << momentumMaximum << endl;
        }
        scalarField bestStep(nkSize,Zero);
        scalar bestLinearNorm=correctionRhsNorm,bestDirectionalDescent=0;
        bool krylovFinite=true,linearConverged=false;
        solveFrozenRhs
        (
            correctionRhs,correctionRhsNorm,stageKrylovRelTol*correctionRhsNorm,
            continuityCorrection ? "coupledContinuity" : "production",
            bestStep,bestLinearNorm,bestDirectionalDescent,krylovFinite,linearConverged
        );
        if (continuityCorrection && krylovFinite && bestLinearNorm<correctionRhsNorm)
        {
            scalarField action(nkSize,Zero);
            krylovFinite=nkApplyJacobian(bestStep,action);
            // Resolve the signed slope in the two physical residual blocks.
            // Positive values are descent for the corresponding squared norm.
            const auto splitDescent = [&](const scalarField& product)
            {
                const scalar continuity=nkDot(continuityRhs,product)/rhsNorm;
                const scalar total=nkDot(rhs,product)/rhsNorm;
                return Pair<scalar>(continuity,total-continuity);
            };
            if (krylovFinite)
            {
                Pair<scalar> descent=splitDescent(action);
                if (stageAuditContinuityCorrection)
                    Info<< "IHBContinuityDescentAudit continuity " << descent.first()
                        << " momentum " << descent.second() << endl;

                if (stageBalanceContinuityMomentum && descent.second()<0
                 && !acceptedMomentumDirection.empty())
                {
                    scalarField momentumAction(nkSize,Zero);
                    if (nkApplyJacobian(acceptedMomentumDirection,momentumAction))
                    {
                        const Pair<scalar> referenceDescent=splitDescent(momentumAction);
                        if (referenceDescent.second()>VSMALL)
                        {
                            // Add the accepted Newton direction only as far as
                            // needed to offset the adverse momentum slope.
                            const scalar factor=-descent.second()/referenceDescent.second();
                            if (std::isfinite(factor))
                            {
                                const scalarField balanced
                                    (bestStep+factor*acceptedMomentumDirection);
                                if (nkApplyJacobian(balanced,momentumAction))
                                {
                                    const Pair<scalar> adjusted=splitDescent(momentumAction);
                                    const bool useBalanced=std::isfinite(adjusted.first())
                                        && std::isfinite(adjusted.second())
                                        && adjusted.first()>0
                                        && adjusted.first()+adjusted.second()>0
                                        && adjusted.second()>descent.second();
                                    if (useBalanced)
                                    {
                                        bestStep=balanced;
                                        descent=adjusted;
                                    }
                                    Info<< "IHBContinuityMomentumBalance factor " << factor
                                        << " referenceContinuity " << referenceDescent.first()
                                        << " referenceMomentum " << referenceDescent.second()
                                        << " continuity " << adjusted.first()
                                        << " momentum " << adjusted.second()
                                        << " used " << useBalanced << endl;
                                }
                            }
                        }
                    }
                }
                bestDirectionalDescent=descent.first()+descent.second();
            }
            if (!std::isfinite(bestDirectionalDescent) || bestDirectionalDescent<=0)
                krylovFinite=false;
            Info<< "IHBStageCoupledContinuity rhs " << correctionRhsNorm
                << " trueFdRelative " << bestLinearNorm/correctionRhsNorm
                << " fullStageDescent " << bestDirectionalDescent << endl;
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
        if (!krylovFinite || bestLinearNorm>=correctionRhsNorm)
        {
            Info<< "IHBStageNewtonRejected newton " << newton
                << " reason invalidOrNonreducingLinearStep"
                << " finite " << krylovFinite << " fdEvaluations " << fdEvaluations
                << " retainedScaledRms " << stageCurrentMerit << endl;
            if (!stageCoupledContinuityCorrection || continuityCorrection) break;
            newton=stageNewtonSteps;
            continue;
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
                << " trueFdRelative " << bestLinearNorm/correctionRhsNorm << endl;
            if (stageEvaluationFinite && trialMerit<=target)
            {
                if (audit)
                {
                    const vector deltaGradient=hotspotGradient()-baseGradient;
                    const scalar beforeDpMax=gMax(mag(nkBasePressure[hotPhase]));
                    const scalar afterDpMax=gMax
                    (
                        mag(scalarSolIncrement[hotPhase].primitiveField())
                    );
                    if (localHotspot)
                        Pout<< "IHBHotspot eta " << eta << " newton " << newton
                            << " beta " << beta << " phase " << hotPhase
                            << " rank " << hotRank << " cell " << hotCell
                            << " centre " << subTimeMeshes[hotPhase].C()[hotCell]
                            << " deltaGradP " << deltaGradient
                            << " seedDpMax " << beforeDpMax
                            << " acceptedDpMax " << afterDpMax << endl;
                }
                accepted=true;
                stageCurrentMerit=trialMerit;
                if (stageBalanceContinuityMomentum && !continuityCorrection)
                    acceptedMomentumDirection=beta*bestStep;
                Info<< "IHBStageNewtonAccepted newton " << newton
                    << " beta " << beta << " scaledRms " << stageCurrentMerit
                    << " previousScaledRms " << nkOldMerit
                    << " fdEvaluations " << fdEvaluations
                    << " rhsKind " << (continuityCorrection ? "coupledContinuity" : "production")
                    << endl;
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
            if (!stageCoupledContinuityCorrection || continuityCorrection) break;
            newton=stageNewtonSteps;
            continue;
        }
        if (continuityCorrection) break;
    }
    if (matchPhysicalDiagonal)
        forAll(subTimeMeshes,K)
            eqSystem.dVByV(K,K).diag()=savedVelocityDiagonal[K];
}
