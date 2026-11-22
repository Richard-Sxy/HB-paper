// Transaction wrapper around the existing epsilon-then-k proposal. Included
// BEFORE building that proposal; phaseAcceptTurbulenceProposal() is called once
// afterwards. The physical k/epsilon equations and their fixed point are unchanged.
PtrList<HBKEpsilon::PhaseImplicitState> phaseTransactionBase(nO);
// Cache the SAME physical residual arrays used for the RMS statistics. Joint
// Newton packs these immediately after measurement, avoiding a second native
// matrix assembly (and repeated wall manipulation) for each FD evaluation.
PtrList<scalarField> phaseNewtonKResidual(nO),phaseNewtonEpsilonResidual(nO);
if (phaseTurbulenceNewton)
    forAll(subTimeMeshes,K)
    {
        phaseNewtonKResidual.set(K,new scalarField(subTimeMeshes[K].nCells(),Zero));
        phaseNewtonEpsilonResidual.set(K,new scalarField(subTimeMeshes[K].nCells(),Zero));
    }
List<boolList> phaseTransactionEpsilonWall(nO);
scalar phaseBaseKRms=0,phaseBaseEpsilonRms=0;
scalar phaseKScale=1,phaseEpsilonScale=1,phaseBaseTurbulenceMerit=0;
bool phaseProposalValid=true;

const auto phaseMeasureTurbulence = [&](scalar& kRms,scalar& epsRms)
{
    bool finite=true;
    forAll(subTimeMeshes,K)
    {
        const HBPhaseTime phaseClock(runTime,snapshots[K]+hbMotionDeltaT,
            hbMotionDeltaT,subTimeMeshes[K],meshPhiHB[K]);
        // k BC -> epsilon wall -> nut -> G and wall production.
        const bool phaseFinite=hbTurbulence[K].preparePhaseImplicitResidual();
        finite=finite && phaseFinite;
    }
    if (!finite) { kRms=GREAT; epsRms=GREAT; return false; }
    // Never carry a candidate's old explicit temporal source into a trial.
    forAll(subTimeMeshes,K)
        hbTurbulence[K].freezeTime(HB,hbTurbulence,K,turbulenceHBDamping);
    scalar kSquare=0,epsSquare=0,volume=0;
    forAll(subTimeMeshes,K)
    {
        const HBPhaseTime phaseClock(runTime,snapshots[K]+hbMotionDeltaT,
            hbMotionDeltaT,subTimeMeshes[K],meshPhiHB[K]);
        const bool phaseFinite=hbTurbulence[K].phaseImplicitResidualStatistics
        (
            kSquare,epsSquare,volume,
            phaseTurbulenceNewton ? &phaseNewtonKResidual[K] : nullptr,
            phaseTurbulenceNewton ? &phaseNewtonEpsilonResidual[K] : nullptr
        );
        finite=finite && phaseFinite;
    }
    reduce(finite,andOp<bool>());
    reduce(kSquare,sumOp<scalar>()); reduce(epsSquare,sumOp<scalar>());
    reduce(volume,sumOp<scalar>());
    finite=finite && std::isfinite(kSquare) && std::isfinite(epsSquare)
        && std::isfinite(volume) && volume>0 && kSquare>=0 && epsSquare>=0;
    if (!finite) { kRms=GREAT; epsRms=GREAT; return false; }
    kRms=Foam::sqrt(kSquare/volume); epsRms=Foam::sqrt(epsSquare/volume);
    return std::isfinite(kRms) && std::isfinite(epsRms);
};
const auto phaseTurbulenceMerit = [&](const scalar kRms,const scalar epsRms)
{
    return sqr(kRms/phaseKScale)+sqr(epsRms/phaseEpsilonScale);
};
const auto phaseRestoreTurbulenceBase = [&]()
{
    forAll(subTimeMeshes,K)
        hbTurbulence[K].restorePhaseImplicitState(phaseTransactionBase[K]);
};
const auto phaseTurbulenceStateDifference = [&]()
{
    scalar difference=0;
    forAll(subTimeMeshes,K)
        difference=max(difference,
            hbTurbulence[K].phaseImplicitStateDifference(phaseTransactionBase[K]));
    reduce(difference,maxOp<scalar>());
    return difference;
};

if (phaseTurbulenceLineSearch)
{
    PtrList<HBKEpsilon::PhaseImplicitState> entryState(nO);
    forAll(subTimeMeshes,K)
    {
        entryState.set(K,new HBKEpsilon::PhaseImplicitState(hbTurbulence[K]));
        phaseTransactionEpsilonWall[K]=hbTurbulence[K].phaseImplicitConstraintMask(true);
    }
    if (!phaseMeasureTurbulence(phaseBaseKRms,phaseBaseEpsilonRms))
        FatalErrorInFunction << "Invalid phaseImplicit canonical baseline residual"
            << exit(FatalError);
    scalar canonicalChange=0,canonicalFieldChange=0;
    forAll(subTimeMeshes,K)
    {
        canonicalChange=max(canonicalChange,
            hbTurbulence[K].phaseImplicitStateDifference(entryState[K]));
        canonicalFieldChange=max(canonicalFieldChange,
            hbTurbulence[K].phaseImplicitStateDifference(entryState[K],false));
        phaseTransactionBase.set(K,new HBKEpsilon::PhaseImplicitState(hbTurbulence[K]));
    }
    reduce(canonicalChange,maxOp<scalar>());
    reduce(canonicalFieldChange,maxOp<scalar>());
    phaseKScale=max(phaseBaseKRms,kTolerance);
    phaseEpsilonScale=max(phaseBaseEpsilonRms,epsilonTolerance);
    phaseBaseTurbulenceMerit=phaseTurbulenceMerit(phaseBaseKRms,phaseBaseEpsilonRms);
    Info<< "IHBPhaseTurbulenceBase kRms " << phaseBaseKRms
        << " epsilonRms " << phaseBaseEpsilonRms
        << " kScale " << phaseKScale << " epsilonScale " << phaseEpsilonScale
        << " merit " << phaseBaseTurbulenceMerit
        << " canonicalFieldChange " << canonicalFieldChange
        << " canonicalNumericStateChangeIncludingCaches " << canonicalChange << endl;
    // Restoring clones clears boundary manipulation flags from measurement,
    // so the next newly assembled proposal receives native wall constraints.
    phaseRestoreTurbulenceBase();
}

const auto phaseAcceptTurbulenceProposal = [&]()
{
    if (!phaseTurbulenceLineSearch) return;
    PtrList<scalarField> proposedK(nO),proposedEpsilon(nO);
    forAll(subTimeMeshes,K)
    {
        proposedK.set(K,new scalarField(hbTurbulence[K].kValues()));
        proposedEpsilon.set(K,new scalarField(hbTurbulence[K].epsilonValues()));
    }
    const auto setTrial = [&](const scalar lambda)
    {
        phaseRestoreTurbulenceBase();
        forAll(subTimeMeshes,K)
        {
            scalarField& k=hbTurbulence[K].phaseImplicitField(false).primitiveFieldRef();
            scalarField& eps=hbTurbulence[K].phaseImplicitField(true).primitiveFieldRef();
            const scalarField& k0=phaseTransactionBase[K].k.primitiveField();
            const scalarField& e0=phaseTransactionBase[K].epsilon.primitiveField();
            // Preserve an EXACT zero step; do not multiply an invalid proposal
            // by zero. Wall epsilon values will be reconstructed from this k.
            if (lambda!=0)
            {
                k=k0+lambda*(proposedK[K]-k0);
                eps=e0+lambda*(proposedEpsilon[K]-e0);
            }
        }
    };
    struct TrialStatistics
    {
        scalar kRelative=0,epsilonRelative=0,minEpsilon=GREAT,maxNut=0;
        scalar worstRelative=0;
        label rank=-1,phase=-1,cell=-1,field=-1;
        bool finite=true;
    };
    const auto trialStatistics = [&]()
    {
        TrialStatistics result;
        scalar localWorst=-1;
        label localPhase=-1,localCell=-1,localField=-1;
        forAll(subTimeMeshes,K)
        {
            const scalarField& k=hbTurbulence[K].kValues();
            const scalarField& eps=hbTurbulence[K].epsilonValues();
            const scalarField& nut=hbTurbulence[K].phaseImplicitNut().primitiveField();
            forAll(k,celli)
            {
                result.finite=result.finite && std::isfinite(k[celli])
                    && std::isfinite(eps[celli]) && std::isfinite(nut[celli]);
                result.minEpsilon=min(result.minEpsilon,eps[celli]);
                result.maxNut=max(result.maxNut,nut[celli]);
                for (label field=0;field<2;++field)
                {
                    // All k unknowns remain free; prescribed epsilon wall
                    // rows are dependent values, not line-search coordinates.
                    if (field==1 && phaseTransactionEpsilonWall[K][celli]) continue;
                    const scalar base=field==0
                        ? phaseTransactionBase[K].k[celli]
                        : phaseTransactionBase[K].epsilon[celli];
                    const scalar value=field==0 ? k[celli] : eps[celli];
                    const scalar relative=mag(value-base)/max
                        (mag(base),hbTurbulence[K].phaseImplicitMinimum(field==1));
                    result.finite=result.finite && std::isfinite(relative);
                    if (field==0) result.kRelative=max(result.kRelative,relative);
                    else result.epsilonRelative=max(result.epsilonRelative,relative);
                    if (relative>localWorst)
                    {
                        localWorst=relative; localPhase=K; localCell=celli; localField=field;
                    }
                }
            }
        }
        reduce(result.finite,andOp<bool>());
        reduce(result.kRelative,maxOp<scalar>());
        reduce(result.epsilonRelative,maxOp<scalar>());
        reduce(result.minEpsilon,minOp<scalar>());
        reduce(result.maxNut,maxOp<scalar>());
        result.worstRelative=max(result.kRelative,result.epsilonRelative);
        label winner=localWorst==result.worstRelative
            ? Pstream::myProcNo() : Pstream::nProcs();
        reduce(winner,minOp<label>());
        scalarList location(3,Zero);
        if (Pstream::myProcNo()==winner)
        {
            location[0]=localPhase; location[1]=localCell; location[2]=localField;
        }
        Pstream::listCombineGather(location,plusEqOp<scalar>());
        Pstream::listCombineScatter(location);
        result.rank=winner; result.phase=label(location[0]);
        result.cell=label(location[1]); result.field=label(location[2]);
        return result;
    };
    const auto logStatistics = [&](const TrialStatistics& values)
    {
        Info<< " kRelativeChange " << values.kRelative
            << " epsilonRelativeChange " << values.epsilonRelative
            << " minEpsilon " << values.minEpsilon << " maxNut " << values.maxNut
            << " limitingRank " << values.rank << " limitingPhase " << values.phase
            << " limitingCell " << values.cell
            << " limitingField " << (values.field==0 ? "k" : "epsilon");
    };
    const auto reject = [&](const word& reason)
    {
        phaseRestoreTurbulenceBase();
        Info<< "IHBPhaseTurbulenceTransaction accepted 0 reason " << reason
            << " lambda 0 kRms " << phaseBaseKRms << " epsilonRms " << phaseBaseEpsilonRms
            << " merit " << phaseBaseTurbulenceMerit
            << " restoredStateDifference " << phaseTurbulenceStateDifference();
        const TrialStatistics restored=trialStatistics();
        logStatistics(restored);
        Info<< endl;
    };

    if (!phaseProposalValid)
    {
        reject("invalidProposal");
        return;
    }

    // Measure the original proposal for comparison, even if the step guard
    // will shorten it. It is not an accepted update merely because it is finite.
    setTrial(1);
    scalar candidateKRms=GREAT,candidateEpsilonRms=GREAT;
    const bool candidateFinite=phaseMeasureTurbulence(candidateKRms,candidateEpsilonRms);
    const scalar candidateMerit=phaseTurbulenceMerit(candidateKRms,candidateEpsilonRms);
    const TrialStatistics candidateStats=trialStatistics();
    Info<< "IHBPhaseTurbulenceCandidate accepted 0 lambda 1 finite " << candidateFinite
        << " kRms " << candidateKRms << " epsilonRms " << candidateEpsilonRms
        << " merit " << candidateMerit;
    logStatistics(candidateStats); Info<< endl;

    // Test restoration FROM the actual proposal, not only an untouched base.
    setTrial(0);
    scalar zeroKRms=GREAT,zeroEpsilonRms=GREAT;
    const bool zeroFinite=phaseMeasureTurbulence(zeroKRms,zeroEpsilonRms);
    const scalar zeroMerit=phaseTurbulenceMerit(zeroKRms,zeroEpsilonRms);
    const scalar zeroResidualError=max(mag(zeroKRms-phaseBaseKRms)/phaseKScale,
        mag(zeroEpsilonRms-phaseBaseEpsilonRms)/phaseEpsilonScale);
    const scalar zeroStateError=phaseTurbulenceStateDifference();
    const scalar zeroMeritError=mag(zeroMerit-phaseBaseTurbulenceMerit);
    Info<< "IHBPhaseTurbulenceZeroStep finite " << zeroFinite
        << " kRms " << zeroKRms << " epsilonRms " << zeroEpsilonRms
        << " scaledResidualError " << zeroResidualError
        << " numericStateErrorIncludingCaches " << zeroStateError << endl;
    phaseRestoreTurbulenceBase();
    if (!zeroFinite || !std::isfinite(zeroResidualError) || zeroResidualError>1e-10
     || !std::isfinite(zeroStateError) || zeroStateError>1e-10)
    {
        reject("zeroStepMismatch");
        return;
    }
    if (!candidateStats.finite || !std::isfinite(candidateStats.worstRelative))
    {
        reject("nonfiniteProposalDirection");
        return;
    }
    if (candidateStats.worstRelative<=VSMALL)
    {
        reject("zeroProposalDirection");
        return;
    }
    const scalar startLambda=min(scalar(1),
        phaseTurbulenceRelativeChangeLimit/candidateStats.worstRelative);
    Info<< "IHBPhaseTurbulenceStepGuard lambda " << startLambda
        << " relativeChangeLimit " << phaseTurbulenceRelativeChangeLimit;
    logStatistics(candidateStats); Info<< endl;
    if (startLambda<phaseTurbulenceMinLambda)
    {
        reject("relativeChangeGuardBelowMinLambda");
        return;
    }

    // A genuine directional Armijo test for the SQUARED joint physical merit.
    // An uncertain/positive slope is rejected, never replaced by a fake descent.
    const scalar probeLambda=min(scalar(1e-4),scalar(1e-3)*startLambda);
    setTrial(probeLambda);
    scalar probeKRms=GREAT,probeEpsilonRms=GREAT;
    const bool probeFinite=phaseMeasureTurbulence(probeKRms,probeEpsilonRms);
    const scalar probeMerit=phaseTurbulenceMerit(probeKRms,probeEpsilonRms);
    const scalar meritChange=probeMerit-phaseBaseTurbulenceMerit;
    const scalar noiseAllowance=10*zeroMeritError
        +1e-10*max(scalar(1),phaseBaseTurbulenceMerit);
    const scalar slope=meritChange/probeLambda;
    Info<< "IHBPhaseTurbulenceSlopeProbe lambda " << probeLambda
        << " finite " << probeFinite << " kRms " << probeKRms
        << " epsilonRms " << probeEpsilonRms << " merit " << probeMerit
        << " meritChange " << meritChange << " noiseAllowance " << noiseAllowance
        << " directionalSlope " << slope << endl;
    if (!probeFinite || !std::isfinite(slope))
    {
        reject("nonfiniteSlopeProbe");
        return;
    }
    if (meritChange>=-noiseAllowance)
    {
        reject(meritChange>noiseAllowance ? "nonDescentDirection" : "uncertainDirectionalSlope");
        return;
    }

    scalar lambda=startLambda;
    for (label trial=0;trial<=phaseTurbulenceMaxBacktracks
         && lambda>=phaseTurbulenceMinLambda;++trial,lambda*=0.5)
    {
        setTrial(lambda);
        scalar trialKRms=GREAT,trialEpsilonRms=GREAT;
        const bool trialFinite=phaseMeasureTurbulence(trialKRms,trialEpsilonRms);
        const scalar trialMerit=phaseTurbulenceMerit(trialKRms,trialEpsilonRms);
        const TrialStatistics stats=trialStatistics();
        const scalar target=phaseBaseTurbulenceMerit+phaseTurbulenceArmijo*lambda*slope;
        const bool accepted=trialFinite && stats.finite && std::isfinite(trialMerit)
            && stats.worstRelative<=phaseTurbulenceRelativeChangeLimit*(1+1e-10)
            && trialMerit<=target;
        Info<< "IHBPhaseTurbulenceTrial trial " << trial << " lambda " << lambda
            << " finite " << trialFinite << " kRms " << trialKRms
            << " epsilonRms " << trialEpsilonRms << " merit " << trialMerit
            << " armijoTarget " << target << " accepted " << accepted;
        logStatistics(stats); Info<< endl;
        if (accepted)
        {
            // Fields, wall state, nut/G and frozen HB sources all belong to
            // this evaluated trial. No stale proposal cache is restored here.
            Info<< "IHBPhaseTurbulenceTransaction accepted 1 lambda " << lambda
                << " kRms " << trialKRms << " epsilonRms " << trialEpsilonRms
                << " merit " << trialMerit << " previousMerit " << phaseBaseTurbulenceMerit
                << " proposalLinearConvergenceIsNotAcceptedStepConvergence 1" << endl;
            return;
        }
    }
    reject("jointResidualBacktracking");
};
