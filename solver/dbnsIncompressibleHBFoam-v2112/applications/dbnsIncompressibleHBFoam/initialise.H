    Info << "Incompressible Harmonic Balance: absolute Cartesian U, mesh-relative phi" << nl << endl;

    PtrList<pseudotimeControl> solnControl(nO);

    forAll(solnControl,i)
    {
    	solnControl.set
    	(
    		i,
    	    new pseudotimeControl
    	    (
    			subTimeMeshes[i],
				true,
				1, // one scalar unknown: kinematic pressure p
				1
    	    )
    	);

        solnControl[i].setCorr(runTime.startTimeIndex());
    }

    Info << "Convective/viscous local pseudo-time selected" << nl << endl;

    #include "createFields.H"

	forAll(subTimeMeshes,i)
	{
	    // Read or initialise pseudo Co number
	    // NOTE: It is not necessarily read before every outer iteration
	    // (see resetPseudo)
		scalar totCells = subTimeMeshes[i].globalData().nTotalCells();

		Info<<"Sub-time mesh "<<i<< nl
			<< "Initial pseudo Courant No: "
			<< "Min: "
			<< min(pseudoCoField[i]).value()
			<< " Mean: "
			<< sum(pseudoCoField[i]).value()/totCells
			<< " Max: "
			<< max(pseudoCoField[i]).value() << nl << endl;
	}


// Physical velocity-pressure stopping criteria, in their own dimensions.
const dictionary& physicalControls = solnControl[0].dict();
const scalar momentumTolerance = physicalControls.lookupOrDefault<scalar>
    ("momentumTolerance", 1e-8); // m/s^2
const scalar continuityTolerance = physicalControls.lookupOrDefault<scalar>
    ("continuityTolerance", 1e-9); // 1/s
const scalar referenceTolerance = physicalControls.lookupOrDefault<scalar>
    ("pressureReferenceTolerance", 1e-9); // m^2/s^2
const scalar fluxTolerance = physicalControls.lookupOrDefault<scalar>
    ("fluxConsistencyTolerance", 1e-10); // |deltaPhi|/|Sf|, m/s
for (const scalar tol : {momentumTolerance, continuityTolerance,
                         referenceTolerance, fluxTolerance})
{
    if (!std::isfinite(tol) || tol <= 0)
        FatalErrorInFunction << "Physical tolerances must be finite and positive"
            << exit(FatalError);
}
scalar momentumMax=GREAT, continuityMax=GREAT;
scalar referenceMax=GREAT, fluxConsistencyMax=GREAT;
bool equationsConverged=false, physicalConverged=false;

const scalar kTolerance=physicalControls.lookupOrDefault<scalar>("kTolerance",1e-8);
const scalar epsilonTolerance=physicalControls.lookupOrDefault<scalar>("epsilonTolerance",1e-8);
const scalar omegaTolerance=physicalControls.lookupOrDefault<scalar>("omegaTolerance",1e-8);
const scalar secondTurbulenceTolerance=hbSST ? omegaTolerance : epsilonTolerance;
const scalar turbulenceHBDamping=physicalControls.lookupOrDefault<scalar>("turbulenceHBDamping",1.0);
const bool turbulenceResidualAudit=physicalControls.lookupOrDefault<bool>
    ("turbulenceResidualAudit",false);
const label nHBTurbulenceCorrectors=physicalControls.lookupOrDefault<label>
    ("nHBTurbulenceCorrectors",1);
if (nHBTurbulenceCorrectors<1)
    FatalErrorInFunction << "nHBTurbulenceCorrectors must be a positive integer"
        << exit(FatalError);
for (const scalar value : {kTolerance,epsilonTolerance,omegaTolerance,turbulenceHBDamping})
    if (!std::isfinite(value) || value<=0)
        FatalErrorInFunction << "HB turbulence tolerances/damping must be finite and positive"
            << exit(FatalError);
scalar kResidualMax=GREAT,epsilonResidualMax=GREAT;
bool turbulenceConverged=!hbRAS;
if (turbulenceHBDamping<1)
    FatalErrorInFunction << "turbulenceHBDamping must be >= 1 for the HB stages"
        << exit(FatalError);

// Optional all-phase implicit scalar turbulence solve. Field increment damping
// is separate from relaxationFactors; this mode never calls equation.relax().
const word hbTurbulenceMethod=physicalControls.lookupOrDefault<word>
    ("hbTurbulenceMethod","ssprk");
const scalar phaseTurbulenceDamping=physicalControls.lookupOrDefault<scalar>
    ("phaseTurbulenceDamping",0.5);
const scalar phaseTurbulenceMinDamping=physicalControls.lookupOrDefault<scalar>
    ("phaseTurbulenceMinDamping",1e-6);
const label phaseTurbulenceDirections=physicalControls.lookupOrDefault<label>
    ("phaseTurbulenceDirections",30);
const label phaseTurbulenceCycles=physicalControls.lookupOrDefault<label>
    ("phaseTurbulenceCycles",5);
const scalar phaseTurbulenceRelTol=physicalControls.lookupOrDefault<scalar>
    ("phaseTurbulenceRelTol",0.01);
// Absolute tolerance is the global volume-weighted RMS of residual / cell V.
const scalar phaseTurbulenceAbsTol=physicalControls.lookupOrDefault<scalar>
    ("phaseTurbulenceAbsTol",1e-12);
const word phaseTurbulencePreconditioner=physicalControls.lookupOrDefault<word>
    ("phaseTurbulencePreconditioner","ihbBlockJacobi");
if (hbTurbulenceMethod!="ssprk" && hbTurbulenceMethod!="phaseImplicit")
    FatalErrorInFunction << "hbTurbulenceMethod must be ssprk or phaseImplicit"
        << exit(FatalError);
if (hbSST && hbTurbulenceMethod!="ssprk")
    FatalErrorInFunction
        << "SST-HB currently uses hbTurbulenceMethod ssprk; phaseImplicit "
        << "is the experimental kEpsilon-only path" << exit(FatalError);
if (!std::isfinite(phaseTurbulenceDamping)
 || phaseTurbulenceDamping<=0 || phaseTurbulenceDamping>1
 || !std::isfinite(phaseTurbulenceMinDamping)
 || phaseTurbulenceMinDamping<=0 || phaseTurbulenceMinDamping>phaseTurbulenceDamping
 || phaseTurbulenceDirections<1 || phaseTurbulenceCycles<1
 || !std::isfinite(phaseTurbulenceRelTol)
 || phaseTurbulenceRelTol<=0 || phaseTurbulenceRelTol>=1
 || !std::isfinite(phaseTurbulenceAbsTol) || phaseTurbulenceAbsTol<=0)
    FatalErrorInFunction << "Invalid phaseImplicit turbulence controls"
        << exit(FatalError);
if (phaseTurbulencePreconditioner!="ihbBlockJacobi"
 && phaseTurbulencePreconditioner!="LUSGS")
    FatalErrorInFunction
        << "phaseTurbulencePreconditioner must be ihbBlockJacobi or LUSGS"
        << exit(FatalError);

// Joint nonlinear acceptance for the experimental phaseImplicit proposal.
// Relative-change protection only shortens the step; it never clips fields.
const bool phaseTurbulenceLineSearch=physicalControls.lookupOrDefault<bool>
    ("phaseTurbulenceLineSearch",true);
const bool phaseTurbulencePseudoTime=physicalControls.lookupOrDefault<bool>
    ("phaseTurbulencePseudoTime",false);
const bool phaseTurbulenceNewton=physicalControls.lookupOrDefault<bool>
    ("phaseTurbulenceNewton",false);
if (phaseTurbulenceNewton
 && (hbTurbulenceMethod!="phaseImplicit" || !phaseTurbulenceLineSearch))
    FatalErrorInFunction << "phaseTurbulenceNewton requires phaseImplicit and its line search"
        << exit(FatalError);
if (phaseTurbulencePseudoTime
 && (hbTurbulenceMethod!="phaseImplicit" || phaseTurbulenceNewton))
    FatalErrorInFunction << "phaseTurbulencePseudoTime applies only to phaseImplicit Picard; "
        << "joint Newton uses the physical Jacobian" << exit(FatalError);
const scalar phaseTurbulenceRelativeChangeLimit=physicalControls.lookupOrDefault<scalar>
    ("phaseTurbulenceRelativeChangeLimit",0.2);
const scalar phaseTurbulenceArmijo=physicalControls.lookupOrDefault<scalar>
    ("phaseTurbulenceArmijo",1e-4);
const label phaseTurbulenceMaxBacktracks=physicalControls.lookupOrDefault<label>
    ("phaseTurbulenceMaxBacktracks",12);
const scalar phaseTurbulenceMinLambda=physicalControls.lookupOrDefault<scalar>
    ("phaseTurbulenceMinLambda",1e-8);
if (!std::isfinite(phaseTurbulenceRelativeChangeLimit)
 || phaseTurbulenceRelativeChangeLimit<=0 || phaseTurbulenceRelativeChangeLimit>=1
 || !std::isfinite(phaseTurbulenceArmijo)
 || phaseTurbulenceArmijo<=0 || phaseTurbulenceArmijo>=0.5
 || phaseTurbulenceMaxBacktracks<0 || phaseTurbulenceMaxBacktracks>30
 || !std::isfinite(phaseTurbulenceMinLambda)
 || phaseTurbulenceMinLambda<=0 || phaseTurbulenceMinLambda>=1)
    FatalErrorInFunction << "Invalid phaseImplicit joint line-search controls"
        << exit(FatalError);

// Optional globalisation of the complete pressure/velocity/face-flux state.
// The merit uses fixed physical scales, never different residual dimensions.
const bool flowLineSearch=physicalControls.lookupOrDefault<bool>("flowLineSearch",false);
// Closed domains: the continuity rows sum to the net boundary flux, which a
// non-conservative AMI makes nonzero. The gauge cell then carries that global
// imbalance; its row is the gauge only and its continuity is reported, not
// required to vanish. false restores the earlier "continuity + gauge" row.
const bool gaugeReplacesContinuity=physicalControls.lookupOrDefault<bool>
    ("gaugeReplacesContinuity",true);
// Pseudo-time continuation of the independent face-closure equation. The
// physical residual/target remain unchanged; only the nonlinear stage changes.
const bool flowFacePseudoTime=physicalControls.lookupOrDefault<bool>("flowFacePseudoTime",false);
const bool stageZeroInitialGuess=physicalControls.lookupOrDefault<bool>("stageZeroInitialGuess",false);
const bool stageAdvectiveMomentum=physicalControls.lookupOrDefault<bool>
    ("stageAdvectiveMomentum",false);
const bool freezeTurbulence=physicalControls.lookupOrDefault<bool>("freezeTurbulence",false);
// Optional fixed-point correction of the omitted interpolated pressure flux.
// This is a diagnostic linearisation refinement, with the same physical flux.
const label pressureFluxCorrectors=
    physicalControls.lookupOrDefault<label>("pressureFluxCorrectors",0);
if (pressureFluxCorrectors<0 || (pressureFluxCorrectors>0 && !flowLineSearch))
    FatalErrorInFunction
        << "pressureFluxCorrectors must be nonnegative and requires flowLineSearch"
        << exit(FatalError);
const scalar flowVelocityScale=physicalControls.lookupOrDefault<scalar>("flowVelocityScale",5.0);
const scalar flowLengthScale=physicalControls.lookupOrDefault<scalar>("flowLengthScale",0.2);
const bool flowArmijoUseSlope=physicalControls.lookupOrDefault<bool>
    ("flowArmijoUseSlope",false);
const scalar flowArmijo=physicalControls.lookupOrDefault<scalar>("flowArmijo",1e-4);
const scalar flowRejectFactor=physicalControls.lookupOrDefault<scalar>("flowRejectFactor",0.5);
const label flowMaxBacktracks=physicalControls.lookupOrDefault<label>("flowMaxBacktracks",8);
if (!std::isfinite(flowVelocityScale) || flowVelocityScale<=0
 || !std::isfinite(flowLengthScale) || flowLengthScale<=0
 || !std::isfinite(flowArmijo) || flowArmijo<=0 || flowArmijo>=1
 || !std::isfinite(flowRejectFactor) || flowRejectFactor<=0 || flowRejectFactor>=1
 || flowMaxBacktracks<0 || flowMaxBacktracks>30)
    FatalErrorInFunction << "Invalid flow line-search controls" << exit(FatalError);
scalar flowMerit=GREAT,flowMomentumRms=GREAT,flowContinuityRms=GREAT,flowFluxRms=GREAT;
bool flowStateFinite=true,flowStepAccepted=true;

// Retain the cell-velocity part of the momentum response to a mass-flux
// correction; the physical residual remains the conservative equation.
const bool massConsistentMomentum=physicalControls.lookupOrDefault<bool>
    ("massConsistentMomentum",false);
if (massConsistentMomentum && !flowLineSearch)
    FatalErrorInFunction << "massConsistentMomentum requires flowLineSearch"
        << exit(FatalError);

// Optional frozen-matrix correction of the COMPLETE nonlinear pseudo-time stage.
// Kept separate from the two partial linearisation experiments above.
const label fullStageCorrectors=
    physicalControls.lookupOrDefault<label>("fullStageCorrectors",0);
const scalar fullStageDamping=
    physicalControls.lookupOrDefault<scalar>("fullStageDamping",1.0);
const label fullStageMaxBacktracks=
    physicalControls.lookupOrDefault<label>("fullStageMaxBacktracks",8);
if (fullStageMaxBacktracks<0 || fullStageMaxBacktracks>30)
    FatalErrorInFunction << "fullStageMaxBacktracks must be in [0,30]"
        << exit(FatalError);
if (fullStageCorrectors<0 || (fullStageCorrectors>0 && !flowLineSearch))
    FatalErrorInFunction
        << "fullStageCorrectors must be nonnegative and requires flowLineSearch"
        << exit(FatalError);
if (!std::isfinite(fullStageDamping) || fullStageDamping<=0 || fullStageDamping>1)
    FatalErrorInFunction << "fullStageDamping must be finite and in (0,1]"
        << exit(FatalError);
if (fullStageCorrectors>0 && (pressureFluxCorrectors>0 || massConsistentMomentum))
    FatalErrorInFunction
        << "fullStageCorrectors is mutually exclusive with pressureFluxCorrectors "
        << "and massConsistentMomentum; the full residual already includes those terms"
        << exit(FatalError);


// Optional matrix-free Newton correction of the same complete frozen-rAU stage.
// Replaces (never follows) the Picard inner loop when enabled.
const label stageNewtonSteps=physicalControls.lookupOrDefault<label>("stageNewtonSteps",0);
// Fraction of the stored outer continuity defect to restore in this stage.
const bool flowProjectCandidateContinuity=physicalControls.lookupOrDefault<bool>
    ("flowProjectCandidateContinuity",false);
if (flowProjectCandidateContinuity && (!stageZeroInitialGuess || flowFacePseudoTime))
    FatalErrorInFunction
        << "Candidate continuity projection requires zero-start stage Newton and physical face response"
        << exit(FatalError);
const bool stageCoupledContinuityCorrection=physicalControls.lookupOrDefault<bool>
    ("stageCoupledContinuityCorrection",false);
if (stageCoupledContinuityCorrection
 && (!stageZeroInitialGuess || flowFacePseudoTime || flowProjectCandidateContinuity))
    FatalErrorInFunction
        << "Coupled continuity correction requires zero-start stage Newton, physical face response and no post-stage projection"
        << exit(FatalError);
const bool stageAuditContinuityCorrection=physicalControls.lookupOrDefault<bool>
    ("stageAuditContinuityCorrection",false);
const bool stageBalanceContinuityMomentum=physicalControls.lookupOrDefault<bool>
    ("stageBalanceContinuityMomentum",false);
if (stageBalanceContinuityMomentum && !stageCoupledContinuityCorrection)
    FatalErrorInFunction << "Momentum-balanced correction requires the coupled continuity stage"
        << exit(FatalError);
const scalar stageContinuityEta=
    physicalControls.lookupOrDefault<scalar>("stageContinuityEta",1.0);
const label stageContinuityStartIndex=
    physicalControls.lookupOrDefault<label>("stageContinuityStartIndex",2);
if (!std::isfinite(stageContinuityEta)
 || stageContinuityEta<=0 || stageContinuityEta>1 || stageContinuityStartIndex<1
 || (stageContinuityEta<1 && (stageNewtonSteps<=0 || flowFacePseudoTime)))
    FatalErrorInFunction << "Invalid stage continuity restoration controls"
        << exit(FatalError);
// Restore this fraction of the outer face-closure defect per stage.
const scalar stageFluxRestoreFraction=
    physicalControls.lookupOrDefault<scalar>("stageFluxRestoreFraction",1.0);
const bool stageLimitFluxRestoreByContinuity=physicalControls.lookupOrDefault<bool>
    ("stageLimitFluxRestoreByContinuity",false);
if (!std::isfinite(stageFluxRestoreFraction)
 || stageFluxRestoreFraction<=0 || stageFluxRestoreFraction>1
 || ((stageFluxRestoreFraction<1 || stageLimitFluxRestoreByContinuity)
  && (stageNewtonSteps<=0 || flowFacePseudoTime)))
    FatalErrorInFunction << "Invalid stage face-closure restoration fraction"
        << exit(FatalError);
// Scale only the continuity rows of the inner-stage nonlinear and Krylov norms.
const scalar stageContinuityWeight=
    physicalControls.lookupOrDefault<scalar>("stageContinuityWeight",1);
if (!std::isfinite(stageContinuityWeight) || stageContinuityWeight<=0)
    FatalErrorInFunction << "stageContinuityWeight must be finite and positive"
        << exit(FatalError);
// The pseudo-time diagonal can be kept solely as preconditioning.
const bool stageMomentumPseudoTime=
    physicalControls.lookupOrDefault<bool>("stageMomentumPseudoTime",true);
if (!stageMomentumPseudoTime && stageNewtonSteps<=0)
    FatalErrorInFunction << "Physical stage momentum requires stage Newton"
        << exit(FatalError);
const label stageKrylovDirections=physicalControls.lookupOrDefault<label>("stageKrylovDirections",30);
const label stageKrylovCycles=physicalControls.lookupOrDefault<label>("stageKrylovCycles",3);
const scalar stageNewtonRelTol=physicalControls.lookupOrDefault<scalar>("stageNewtonRelTol",0.1);
const scalar stageKrylovRelTol=physicalControls.lookupOrDefault<scalar>("stageKrylovRelTol",0.1);
const scalar stageFiniteDifferenceEpsilon=
    physicalControls.lookupOrDefault<scalar>("stageFiniteDifferenceEpsilon",1e-7);
if (stageNewtonSteps<0 || stageKrylovDirections<1 || stageKrylovCycles<1
 || !std::isfinite(stageNewtonRelTol) || stageNewtonRelTol<=0 || stageNewtonRelTol>=1
 || !std::isfinite(stageKrylovRelTol) || stageKrylovRelTol<=0 || stageKrylovRelTol>=1
 || !std::isfinite(stageFiniteDifferenceEpsilon)
 || stageFiniteDifferenceEpsilon<=0 || stageFiniteDifferenceEpsilon>=1)
    FatalErrorInFunction << "Invalid matrix-free stage Newton/Krylov controls"
        << exit(FatalError);
if (stageNewtonSteps>0 && (fullStageCorrectors<=0 || !flowLineSearch))
    FatalErrorInFunction
        << "stageNewtonSteps requires fullStageCorrectors>0 and flowLineSearch"
        << exit(FatalError);
if ((flowFacePseudoTime || stageZeroInitialGuess)
 && (!flowLineSearch || stageNewtonSteps<=0 || pressureFluxCorrectors>0 || massConsistentMomentum))
    FatalErrorInFunction
        << "Face pseudo-time/zero initial guess requires stage Newton and no partial flux correctors"
        << exit(FatalError);
const bool flowRetryWithAlternateStage=
    physicalControls.lookupOrDefault<bool>("flowRetryWithAlternateStage",false);
const bool flowRetryWithShiftedStage=
    physicalControls.lookupOrDefault<bool>("flowRetryWithShiftedStage",false);
if (flowRetryWithShiftedStage && !stageZeroInitialGuess)
    FatalErrorInFunction << "Retrying with a different pseudo-time term requires the zero-start stage Newton path"
        << exit(FatalError);
if ((stageAdvectiveMomentum || flowRetryWithAlternateStage) && !stageZeroInitialGuess)
    FatalErrorInFunction << "An advective momentum stage requires the zero-start stage Newton path"
        << exit(FatalError);

// Pressure Schur is an optional RIGHT preconditioner for the same matrix-free J.
// GAMG may vary its coarse-grid action; stored Z vectors implement FGMRES.
const word stagePreconditioner=physicalControls.lookupOrDefault<word>
    ("stagePreconditioner","blockJacobi");
const label stageSchurVCycles=physicalControls.lookupOrDefault<label>
    ("stageSchurVCycles",2);
const scalar stageSchurPivotTolerance=physicalControls.lookupOrDefault<scalar>
    ("stageSchurPivotTolerance",1e-12);
if (stagePreconditioner!="blockJacobi" && stagePreconditioner!="pressureSchur")
    FatalErrorInFunction << "stagePreconditioner must be blockJacobi or pressureSchur"
        << exit(FatalError);
if (stageSchurVCycles<1 || stageSchurVCycles>8
 || !std::isfinite(stageSchurPivotTolerance)
 || stageSchurPivotTolerance<=0 || stageSchurPivotTolerance>=1)
    FatalErrorInFunction << "Require stageSchurVCycles in [1,8] and "
        << "stageSchurPivotTolerance in (0,1)" << exit(FatalError);

// One internal-face block LU-SGS pair; zero preserves the local D inverse.
const label stageSchurVelocitySweeps=physicalControls.lookupOrDefault<label>
    ("stageSchurVelocitySweeps",0);
const label stageSchurVelocityCorrections=
    physicalControls.lookupOrDefault<label>("stageSchurVelocityCorrections",0);
if (stageSchurVelocityCorrections<0
 || (stageSchurVelocityCorrections>0 && stagePreconditioner!="pressureSchur"))
    FatalErrorInFunction << "Velocity defect corrections require pressureSchur and a nonnegative count"
        << exit(FatalError);
const bool stageSchurCorrectFlux=
    physicalControls.lookupOrDefault<bool>("stageSchurCorrectFlux",false);
if (stageSchurCorrectFlux && stagePreconditioner!="pressureSchur")
    FatalErrorInFunction << "Schur flux correction requires pressureSchur"
        << exit(FatalError);
const bool stageAuditVelocityBlock=physicalControls.lookupOrDefault<bool>
    ("stageAuditVelocityBlock",false);
if (stageAuditVelocityBlock && (stagePreconditioner!="pressureSchur" || stageNewtonSteps<=0))
    FatalErrorInFunction << "Velocity block audit requires pressureSchur stage Newton"
        << exit(FatalError);
const bool stageVelocityAuditMatchFlux=physicalControls.lookupOrDefault<bool>
    ("stageVelocityAuditMatchFlux",false);
const bool stageAuditConvectionFlux=physicalControls.lookupOrDefault<bool>
    ("stageAuditConvectionFlux",false);
const bool stageSchurVelocityConvectionFeedback=
    physicalControls.lookupOrDefault<bool>("stageSchurVelocityConvectionFeedback",false);
const bool stageSchurConvectionFeedback=physicalControls.lookupOrDefault<bool>
    ("stageSchurConvectionFeedback",false);
if (stageSchurConvectionFeedback && stagePreconditioner!="pressureSchur")
    FatalErrorInFunction << "Convection feedback requires pressureSchur"
        << exit(FatalError);
const bool stageSchurMatchMomentumDiagonal=physicalControls.lookupOrDefault<bool>
    ("stageSchurMatchMomentumDiagonal",false);
if (stageSchurMatchMomentumDiagonal && stagePreconditioner!="pressureSchur")
    FatalErrorInFunction << "Matched momentum diagonal requires pressureSchur"
        << exit(FatalError);
const bool stageSchurCheckSweeps=physicalControls.lookupOrDefault<bool>
    ("stageSchurCheckSweeps",false);
if (stageSchurVelocitySweeps<0 || stageSchurVelocitySweeps>1)
    FatalErrorInFunction << "stageSchurVelocitySweeps must be 0 or 1"
        << exit(FatalError);
const scalar stageNewtonArmijo=physicalControls.lookupOrDefault<scalar>
    ("stageNewtonArmijo",1e-4);
if (!std::isfinite(stageNewtonArmijo) || stageNewtonArmijo<=0 || stageNewtonArmijo>=0.5)
    FatalErrorInFunction << "stageNewtonArmijo must be finite and in (0,0.5)"
        << exit(FatalError);
