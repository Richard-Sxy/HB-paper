// Scalar-only exact-matrix RIGHT FGMRES for one frozen all-phase equation.
// Included after phaseMatrix and its independent native-RHS check exist.
// Ahat*x=(M*x)/V, bhat=r/V; dot uses global volume weights. The right PC
// applies P^-1*(V*v), so its output is a physical field increment. All Z
// vectors are retained; no shared GMRES solve or finite differences are used.
{
    const auto ptFinite = [&](const scalarField& values)
    {
        bool finite=true;
        forAll(values,i) finite=finite && std::isfinite(values[i]);
        reduce(finite,andOp<bool>());
        return finite;
    };
    const auto ptDot = [&](const scalarField& a,const scalarField& b)
    {
        scalar value=0;
        forAll(a,i) value+=phaseWeights[i]*a[i]*b[i];
        reduce(value,sumOp<scalar>());
        if (!std::isfinite(value))
            FatalErrorInFunction << "Nonfinite phaseImplicit Krylov inner product for "
                << phaseFieldName << exit(FatalError);
        return value;
    };
    const auto ptNorm = [&](const scalarField& a)
    {
        const scalar square=ptDot(a,a);
        if (square<0)
            FatalErrorInFunction << "Negative phaseImplicit norm square"
                << exit(FatalError);
        return Foam::sqrt(square);
    };
    const auto ptApplyMatrix = [&](const scalarField& input,scalarField& output)
    {
        if (!ptFinite(input))
            FatalErrorInFunction << "Nonfinite phaseImplicit matrix input for "
                << phaseFieldName << exit(FatalError);
        forAll(subTimeMeshes,K)
        {
            scalarField& values=phaseWorkFields[K].primitiveFieldRef();
            forAll(values,celli) values[celli]=input[phaseOffsets[K]+celli];
        }
        // Physical patch fields are calculated and never evaluated: their
        // homogeneous action is already in M.diag(). Coupled patch types use
        // each phase's own mesh/AMI addressing through these scratch fields.
        phaseMatrix.matrixMul(phaseWorkFields,phaseNoVectorFields,
            phaseMatrixAction,phaseNoVectorWork);
        output.setSize(phaseFlatSize);
        forAll(subTimeMeshes,K)
        {
            const scalarField& volumes=subTimeMeshes[K].V();
            forAll(volumes,celli)
                output[phaseOffsets[K]+celli]=phaseMatrixAction[K][celli]/volumes[celli];
        }
        if (!ptFinite(output))
            FatalErrorInFunction << "Nonfinite phaseImplicit matrix action for "
                << phaseFieldName << exit(FatalError);
    };

    if (!ptFinite(phaseRhs))
        FatalErrorInFunction << "Nonfinite phaseImplicit RHS for " << phaseFieldName
            << exit(FatalError);
    const scalar rhsNorm=ptNorm(phaseRhs);
    const scalar linearTarget=max(phaseTurbulenceAbsTol,phaseTurbulenceRelTol*rhsNorm);
    scalar bestNorm=rhsNorm;
    label totalDirections=0,completedCycles=0;
    bool stagnated=false;
    Info<< "IHBPhaseTurbulenceLinearStart field " << phaseFieldName
        << " method FGMRES side right operator exactM residualScaling inverseCellVolume"
        << " preconditioner " << phaseTurbulencePreconditioner
        << " initialRms " << rhsNorm << " targetRms " << linearTarget
        << " directions " << phaseTurbulenceDirections
        << " cycles " << phaseTurbulenceCycles << endl;

    if (rhsNorm>linearTarget)
    {
        // A solver instance supplies the PUBLIC owner API to the runtime PC.
        // Neither owner->solve nor owner->solveDelta is ever called.
        dictionary controls;
        controls.add("nDirections",phaseTurbulenceDirections);
        controls.add("maxIter",phaseTurbulenceCycles);
        controls.add("tolerance",phaseTurbulenceAbsTol);
        controls.add("relTol",phaseTurbulenceRelTol);
        controls.add("preconditioner",phaseTurbulencePreconditioner);
        controls.add("blockJacobiCheckInverse",true);
        dictionary ownerControls;
        ownerControls.add("solver",word("GMRES"));
        ownerControls.add("GMRES",controls);
        autoPtr<coupledMatrix::solver> owner=
            coupledMatrix::solver::New(ownerControls,phaseMatrix);
        autoPtr<coupledMatrix::preconditioner> preconditioner=
            coupledMatrix::preconditioner::New(owner(),controls);
        const auto ptApplyPreconditioner = [&](const scalarField& input,scalarField& output)
        {
            forAll(subTimeMeshes,K)
            {
                const scalarField& volumes=subTimeMeshes[K].V();
                forAll(volumes,celli)
                    phasePcWork[K][celli]=volumes[celli]*input[phaseOffsets[K]+celli];
            }
            preconditioner->precondition(phasePcWork,phaseNoVectorWork);
            output.setSize(phaseFlatSize);
            forAll(subTimeMeshes,K)
                forAll(phasePcWork[K],celli)
                    output[phaseOffsets[K]+celli]=phasePcWork[K][celli];
            if (!ptFinite(output))
                FatalErrorInFunction << "Nonfinite phaseImplicit preconditioner action for "
                    << phaseFieldName << exit(FatalError);
        };

        scalarField totalStep(phaseFlatSize,Zero),residual(phaseRhs);
        for (label cycle=0;cycle<phaseTurbulenceCycles;++cycle)
        {
            const scalar beta=ptNorm(residual);
            if (beta<=linearTarget) break;
            const label m=phaseTurbulenceDirections;
            List<scalarField> basis(m+1),preconditioned(m),hessenberg(m+1);
            for (label row=0;row<=m;++row) hessenberg[row]=scalarField(m,Zero);
            scalarField cosine(m,Zero),sine(m,Zero),g(m+1,Zero);
            basis[0]=residual/beta;
            g[0]=beta;
            label used=0;
            bool happyBreakdown=false,singularDirection=false;
            scalar estimatedNorm=beta;
            for (label column=0;column<m;++column)
            {
                ptApplyPreconditioner(basis[column],preconditioned[column]);
                scalarField work(phaseFlatSize,Zero);
                ptApplyMatrix(preconditioned[column],work);
                const scalar originalNorm=ptNorm(work);
                // Two-pass modified Gram-Schmidt in the same global metric.
                for (label sweep=0;sweep<2;++sweep)
                    for (label row=0;row<=column;++row)
                    {
                        const scalar projection=ptDot(work,basis[row]);
                        hessenberg[row][column]+=projection;
                        work-=projection*basis[row];
                    }
                const scalar nextNorm=ptNorm(work);
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
                if (!std::isfinite(rotationNorm))
                    FatalErrorInFunction << "Nonfinite phaseImplicit Givens rotation"
                        << exit(FatalError);
                if (rotationNorm<=VSMALL)
                {
                    singularDirection=true;
                    break; // No zero division; retain previous usable columns.
                }
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
            totalDirections+=used;
            ++completedCycles;
            if (!used)
            {
                stagnated=true;
                Info<< "IHBPhaseTurbulenceKrylov field " << phaseFieldName
                    << " cycle " << cycle+1 << " directions 0 singularDirection 1"
                    << " trueRelative " << bestNorm/rhsNorm << endl;
                break;
            }
            scalarField coefficients(used,Zero);
            for (label row=used-1;row>=0;--row)
            {
                scalar value=g[row];
                for (label col=row+1;col<used;++col)
                    value-=hessenberg[row][col]*coefficients[col];
                const scalar diagonal=hessenberg[row][row];
                if (!std::isfinite(value) || !std::isfinite(diagonal) || mag(diagonal)<=VSMALL)
                    FatalErrorInFunction << "Invalid phaseImplicit triangular Krylov solve"
                        << exit(FatalError);
                coefficients[row]=value/diagonal;
            }
            for (label column=0;column<used;++column)
                totalStep+=coefficients[column]*preconditioned[column];
            if (!ptFinite(totalStep))
                FatalErrorInFunction << "Nonfinite phaseImplicit Krylov correction"
                    << exit(FatalError);
            scalarField trueAction(phaseFlatSize,Zero);
            ptApplyMatrix(totalStep,trueAction);
            residual=phaseRhs-trueAction;
            const scalar trueNorm=ptNorm(residual);
            if (trueNorm<bestNorm)
            {
                bestNorm=trueNorm;
                phaseIncrement=totalStep;
            }
            Info<< "IHBPhaseTurbulenceKrylov field " << phaseFieldName
                << " cycle " << cycle+1 << " directions " << used
                << " happyBreakdown " << happyBreakdown
                << " singularDirection " << singularDirection
                << " estimatedRelative " << estimatedNorm/rhsNorm
                << " trueRelative " << trueNorm/rhsNorm
                << " targetRelative " << phaseTurbulenceRelTol << endl;
            if (trueNorm<=linearTarget) break;
            if (happyBreakdown || singularDirection)
            {
                stagnated=true;
                break;
            }
            // Restart from the best explicitly verified finite correction.
            if (trueNorm>bestNorm)
            {
                totalStep=phaseIncrement;
                ptApplyMatrix(totalStep,trueAction);
                residual=phaseRhs-trueAction;
            }
        }
    }

    // Native epsilon wall constraints prescribe the CURRENT internal value:
    // their increment is exactly zero. Verify the actual, projected correction
    // independently against M, including every interface, before field damping.
    forAll(subTimeMeshes,K)
        forAll(phaseConstraintMask[K],celli)
            if (phaseConstraintMask[K][celli])
                phaseIncrement[phaseOffsets[K]+celli]=0;
    scalarField finalAction(phaseFlatSize,Zero);
    ptApplyMatrix(phaseIncrement,finalAction);
    const scalarField finalResidual(finalAction-phaseRhs);
    const scalar finalNorm=ptNorm(finalResidual);
    phaseLinearRelative=rhsNorm>VSMALL ? finalNorm/rhsNorm : scalar(0);
    phaseLinearConverged=finalNorm<=linearTarget;
    Info<< "IHBPhaseTurbulenceLinearEnd field " << phaseFieldName
        << " initialRms " << rhsNorm << " finalRms " << finalNorm
        << " trueRelative " << phaseLinearRelative
        << " targetRms " << linearTarget
        << " linearConverged " << phaseLinearConverged
        << " totalDirections " << totalDirections << " completedCycles " << completedCycles
        << " rhsBelowAbsoluteTolerance " << (rhsNorm<=phaseTurbulenceAbsTol)
        << " stagnated " << stagnated << endl;
    forAll(subTimeMeshes,K)
    {
        scalar maximum=0;
        forAll(phaseBase[K],celli)
            maximum=max(maximum,mag(finalResidual[phaseOffsets[K]+celli]));
        reduce(maximum,maxOp<scalar>());
        Info<< "IHBPhaseTurbulenceLinearResidual field " << phaseFieldName
            << " phase " << K << " maxAbsMdeltaMinusRoverV " << maximum << endl;
    }
    if (!phaseLinearConverged)
        Info<< "IHBPhaseTurbulenceLinearTargetNotReached field " << phaseFieldName
            << " retainedFiniteBestCorrection 1 trueRelative " << phaseLinearRelative
            << " targetRms " << linearTarget << " finalRms " << finalNorm << endl;
}
