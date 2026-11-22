scalarField velocityRhs(rhs);
for (label i=0;i<nkSize;i+=4) velocityRhs[i]=0;
const scalar velocityRhsNorm=nkNorm(velocityRhs);
if (velocityRhsNorm>VSMALL)
{
    const VelocityProbe modes[] =
    {
        VelocityProbe::assembled, VelocityProbe::frozenFlux,
        VelocityProbe::frozenFlux, VelocityProbe::coupled
    };
    const char* names[] =
    {
        matchedAudit ? "auditAssembledMatchedFlux" : "auditAssembledVelocity",
        "auditFrozenFluxVelocity", "auditReconstructionVelocity",
        "auditCoupledMomentum"
    };
    for (label mode=0;mode<4;++mode)
    {
        if (mode==2 && !matchedAudit) continue;
        if (mode==3 && matchedAudit) continue;
        velocityProbe=modes[mode];
        applyVelocityReconstructionFeedback=mode==2;
        scalarField auditStep(nkSize,Zero),action(nkSize,Zero);
        scalar bestNorm=velocityRhsNorm,descent=0;
        bool finite=true,converged=false;
        const label before=fdEvaluations;
        solveFrozenRhs
        (
            velocityRhs,velocityRhsNorm,stageKrylovRelTol*velocityRhsNorm,names[mode],
            auditStep,bestNorm,descent,finite,converged
        );
        if (!nkApplyJacobian(auditStep,action)) finite=false;
        const scalarField residual(velocityRhs-action);
        const scalar total=nkNorm(residual);
        Info<< "IHBVelocityBlockAudit operator "
            << names[mode]
            << " initialMomentumRms " << velocityRhsNorm
            << " momentumRelative "
            << Foam::sqrt(nkMomentumDot(residual,residual))/velocityRhsNorm
            << " totalRelative " << total/velocityRhsNorm
            << " finite " << finite
            << " converged " << converged
            << " fdEvaluations " << fdEvaluations-before << endl;
        if (matchedAudit && mode==0)
        {
            // Measure the assembled solution under the true frozen-flux
            // derivative at the very same base. Both actions keep p,q fixed.
            velocityProbe=VelocityProbe::frozenFlux;
            captureVelocityMismatch=true;
            const bool actionValid=nkApplyJacobian(auditStep,action);
            captureVelocityMismatch=false;
            const scalarField trueResidual(velocityRhs-action);
            Info<< "IHBVelocityBlockMatchedTrueResidual relative "
                << nkNorm(trueResidual)/velocityRhsNorm
                << " finite " << actionValid << endl;
            if (actionValid)
            {
                const scalar rhsSquare=sqr(velocityRhsNorm);
                const scalarField sumResidual(trueResidual+residual);
                scalarField reconstruction(nkSize,Zero);
                const char* terms[]={"convection","viscousStress","HB"};
                for (label term=0;term<3;++term)
                {
                    const scalarField& part=velocityBlockMismatch[term];
                    reconstruction+=part;
                    Info<< "IHBVelocityDefectTerm " << terms[term]
                        << " mismatchRelative " << nkNorm(part)/velocityRhsNorm
                        << " residualSquareChange "
                        << -nkDot(sumResidual,part)/rhsSquare << endl;
                }
                const scalarField convectionCorrectedResidual
                (
                    residual-velocityBlockMismatch[0]
                );
                Info<< "IHBVelocityReconstructionOnAssembledStep relative "
                    << nkNorm(convectionCorrectedResidual)/velocityRhsNorm << endl;
                Info<< "IHBVelocityDefectClosure"
                    << " mismatchRelative " << nkNorm(scalarField(residual-trueResidual))/velocityRhsNorm
                    << " reconstructionError " << nkNorm
                    (
                        scalarField(residual-trueResidual-reconstruction)
                    )/velocityRhsNorm
                    << " residualSquareChange "
                    << (nkDot(trueResidual,trueResidual)-nkDot(residual,residual))/rhsSquare
                    << endl;
            }
        }
    }
    reduce(maxFrozenMassAction,maxOp<scalar>());
    applyVelocityReconstructionFeedback=false;
    velocityProbe=VelocityProbe::production;
    nkRestoreBase();
    evaluateFullStage(newton,false);
    const scalar restoredError=nkNorm(scalarField(rhs-nkPackNegativeH()))/rhsNorm;
    Info<< "IHBVelocityBlockAuditRestore relativeError " << restoredError
        << " maxFrozenMassAction " << maxFrozenMassAction << endl;
    if (!stageEvaluationFinite || !std::isfinite(restoredError)
     || restoredError<0 || restoredError>1e-10)
        FatalErrorInFunction << "Velocity-block audit changed the Newton baseline"
            << exit(FatalError);
    fdEvaluations=0;
}
