// All-phase turbulence update: epsilon-then-k Picard or joint physical Newton.
// Included in the outer turbulence-corrector loop; U and phi stay frozen.
// Joint Newton rebuilds nut, production and wall constraints for every trial.
// The native spatial matrix, physical boundary source and coupled interfaces
// are inserted once. Only HB phase coupling is added by this application.
{
    struct phaseTurbulenceFloatingPointGuard
    {
        sigFpe::ignore ignore;
        ~phaseTurbulenceFloatingPointGuard()
        {
            std::feclearexcept(FE_ALL_EXCEPT);
        }
    } phaseArithmeticGuard;

    const fvMesh& phaseReferenceMesh=subTimeMeshes[0];
    bool topologyMatches=true;
    forAll(subTimeMeshes,K)
    {
        hbTurbulence[K].validatePhaseImplicit();
        const fvMesh& meshK=subTimeMeshes[K];
        bool matches=meshK.nCells()==phaseReferenceMesh.nCells()
            && meshK.nFaces()==phaseReferenceMesh.nFaces()
            && meshK.nInternalFaces()==phaseReferenceMesh.nInternalFaces()
            && meshK.boundary().size()==phaseReferenceMesh.boundary().size()
            && meshK.comm()==phaseReferenceMesh.comm();
        if (matches)
        {
            forAll(meshK.faceOwner(),facei)
                matches=matches
                    && meshK.faceOwner()[facei]==phaseReferenceMesh.faceOwner()[facei];
            forAll(meshK.faceNeighbour(),facei)
                matches=matches
                    && meshK.faceNeighbour()[facei]==phaseReferenceMesh.faceNeighbour()[facei];
            forAll(meshK.boundary(),patchi)
            {
                const fvPatch& patch=meshK.boundary()[patchi];
                const fvPatch& reference=phaseReferenceMesh.boundary()[patchi];
                const bool patchMatches=patch.name()==reference.name()
                    && patch.type()==reference.type() && patch.size()==reference.size();
                matches=matches && patchMatches;
                if (patchMatches)
                    forAll(patch.faceCells(),facei)
                        matches=matches && patch.faceCells()[facei]==reference.faceCells()[facei];
            }
        }
        topologyMatches=topologyMatches && matches;
    }
    reduce(topologyMatches,andOp<bool>());
    if (!topologyMatches)
        FatalErrorInFunction
            << "phaseImplicit requires corresponding cell, face and patch topology "
            << "on every phase mesh and MPI rank" << exit(FatalError);

    // Like native kEpsilon, production is prepared before epsilon, while the
    // k destruction coefficient sees the newly updated epsilon. Update nut only
    // after both all-phase solves, preserving this segregated model ordering.
    forAll(subTimeMeshes,K)
    {
        const HBPhaseTime phaseClock(runTime,snapshots[K]+hbMotionDeltaT,
            hbMotionDeltaT,subTimeMeshes[K],meshPhiHB[K]);
        hbTurbulence[K].prepare();
    }

    #include "phaseTurbulenceLineSearch.H"

    if (phaseTurbulenceNewton)
    {
        #include "phaseTurbulenceNewton.H"
    }
    else
    for (const bool phaseSolveEpsilon : {true,false})
    {
        const word phaseFieldName(phaseSolveEpsilon ? "epsilon" : "k");
        List<boolList> phaseConstraintMask(nO);
        PtrList<scalarField> phaseHbDiagonal(nO*nO);
        PtrList<fvScalarMatrix> phaseNativeEquation(nO);
        PtrList<scalarField> phaseBase(nO),phaseExpectedRhs(nO),phaseRhsScale(nO);
        forAll(subTimeMeshes,K)
        {
            const HBPhaseTime phaseClock(runTime,snapshots[K]+hbMotionDeltaT,
                hbMotionDeltaT,subTimeMeshes[K],meshPhiHB[K]);
            phaseConstraintMask[K]=hbTurbulence[K].phaseImplicitConstraintMask(phaseSolveEpsilon);
            forAll(subTimeMeshes,J)
            {
                blockFvMatrix<scalar,scalar> hbBlock(subTimeMeshes[K]);
                hbBlock.diag()=Zero;
                HB.addBlock(hbBlock,K,J);
                phaseHbDiagonal.set(K*nO+J,new scalarField(hbBlock.diag()));
                scalarField& vd=phaseHbDiagonal[K*nO+J];
                forAll(vd,celli)
                    if (phaseConstraintMask[K][celli]) vd[celli]=0;
            }
            tmp<fvScalarMatrix> native=hbTurbulence[K].phaseImplicitEquation(phaseSolveEpsilon);
            native.ref().diag()+=phaseHbDiagonal[K*nO+K];
            if (phaseTurbulencePseudoTime)
            {
                // V/dtau * (q-qOld) regularizes the increment only. Adding
                // the matching source leaves the physical RHS unchanged.
                // Use the SAME local pseudo-time as momentum; HB remains
                // implicit across all phases. Apply before native wall rows.
                const scalarField pseudoDiagonal
                    (subTimeMeshes[K].V()*rPseudoDeltaT[K].primitiveField());
                native.ref().diag()+=pseudoDiagonal;
                native.ref().source()+=pseudoDiagonal
                    *hbTurbulence[K].phaseImplicitField(phaseSolveEpsilon).primitiveField();
            }
            // No relax(), no SSPRK timeRate, no explicit old HB source.
            // Wall setValues replaces the COMPLETE row, including HB self.
            hbTurbulence[K].phaseImplicitConstrain(native.ref(),phaseSolveEpsilon);
            phaseNativeEquation.set(K,native.ptr());
        }
        // Constraints may update q, so take all snapshots after all phases have
        // applied their native constraints, before inserting any native source.
        forAll(subTimeMeshes,K)
            phaseBase.set(K,new scalarField
                (hbTurbulence[K].phaseImplicitField(phaseSolveEpsilon).primitiveField()));
        forAll(subTimeMeshes,K)
        {
            const HBPhaseTime phaseClock(runTime,snapshots[K]+hbMotionDeltaT,
                hbMotionDeltaT,subTimeMeshes[K],meshPhiHB[K]);
            const volScalarField nativeH(phaseNativeEquation[K].H());
            const volScalarField nativeA(phaseNativeEquation[K].A());
            const scalarField& volumes=subTimeMeshes[K].V();
            phaseExpectedRhs.set(K,new scalarField
                (volumes*(nativeH.primitiveField()-nativeA.primitiveField()*phaseBase[K])));
            phaseRhsScale.set(K,new scalarField
                (volumes*(mag(nativeH.primitiveField())
                    +mag(nativeA.primitiveField()*phaseBase[K]))));
            forAll(subTimeMeshes,J)
                if (J!=K)
                {
                    phaseExpectedRhs[K]-=phaseHbDiagonal[K*nO+J]*phaseBase[J];
                    phaseRhsScale[K]+=mag(phaseHbDiagonal[K*nO+J]*phaseBase[J]);
                }
        }

        coupledMatrix phaseMatrix(subTimeMeshes[0],nO,0,true);
        phaseMatrix.setEquationMeshes
            (hbPhaseMeshes(subTimeMeshes),UPtrList<const fvMesh>(0));
        forAll(subTimeMeshes,K)
            forAll(subTimeMeshes,J)
            {
                blockFvMatrix<scalar,scalar>& block=phaseMatrix.dSByS(K,J);
                block.diag()=Zero;
                block.source()=Zero;
                if (K==J)
                {
                    // insertEquation adds physical boundary source ONCE and
                    // copies native coupled interface coefficients. It mutates
                    // native.source(), hence the independent check above.
                    block.insertEquation(phaseNativeEquation[K]);
                }
                else block.diag()=phaseHbDiagonal[K*nO+J];
            }

        labelList phaseOffsets(nO+1,0);
        scalar phaseVolume=0;
        bool phaseStateFinite=true;
        PtrList<volScalarField> phaseWorkFields(nO);
        PtrList<scalarField> phaseMatrixAction(nO),phaseIntegratedRhs(nO),phasePcWork(nO);
        PtrList<volVectorField> phaseNoVectorFields(0);
        PtrList<vectorField> phaseNoVectorWork(0);
        forAll(subTimeMeshes,K)
        {
            const volScalarField& q=hbTurbulence[K].phaseImplicitField(phaseSolveEpsilon);
            const scalar minimum=hbTurbulence[K].phaseImplicitMinimum(phaseSolveEpsilon);
            wordList patchTypes(q.boundaryField().types());
            forAll(patchTypes,patchi)
                if (!q.boundaryField()[patchi].coupled() && patchTypes[patchi]!="empty")
                    patchTypes[patchi]="calculated";
            phaseWorkFields.set(K,new volScalarField
            (
                IOobject("phaseImplicitWork_"+phaseFieldName+Foam::name(K),
                    runTime.timeName(),subTimeMeshes[K],IOobject::NO_READ,
                    IOobject::NO_WRITE,false),
                subTimeMeshes[K],dimensionedScalar("zero",q.dimensions(),0),patchTypes
            ));
            phaseWorkFields[K].primitiveFieldRef()=phaseBase[K];
            phaseMatrixAction.set(K,new scalarField(q.size(),Zero));
            phaseIntegratedRhs.set(K,new scalarField(q.size(),Zero));
            phasePcWork.set(K,new scalarField(q.size(),Zero));
            phaseOffsets[K+1]=phaseOffsets[K]+q.size();
            const scalarField& volumes=subTimeMeshes[K].V();
            phaseStateFinite=phaseStateFinite && std::isfinite(minimum) && minimum>0;
            forAll(volumes,celli)
            {
                phaseVolume+=volumes[celli];
                phaseStateFinite=phaseStateFinite && std::isfinite(volumes[celli])
                    && volumes[celli]>0 && std::isfinite(phaseBase[K][celli])
                    && phaseBase[K][celli]>=minimum;
            }
        }
        reduce(phaseStateFinite,andOp<bool>());
        reduce(phaseVolume,sumOp<scalar>());
        if (!phaseStateFinite || !std::isfinite(phaseVolume) || phaseVolume<=0)
            FatalErrorInFunction << "Invalid phaseImplicit initial " << phaseFieldName
                << " field, model minimum, or cell volume" << exit(FatalError);
        const label phaseFlatSize=phaseOffsets[nO];
        scalarField phaseWeights(phaseFlatSize,Zero),phaseRhs(phaseFlatSize,Zero);
        phaseMatrix.matrixMul(phaseWorkFields,phaseNoVectorFields,
            phaseMatrixAction,phaseNoVectorWork);
        scalar assemblyBackwardError=0,assemblyDifferenceSquare=0,assemblyReferenceSquare=0;
        bool assemblyFinite=true;
        forAll(subTimeMeshes,K)
        {
            scalarField source(subTimeMeshes[K].nCells(),Zero);
            forAll(subTimeMeshes,J) source+=phaseMatrix.dSByS(K,J).source();
            phaseIntegratedRhs[K]=source-phaseMatrixAction[K];
            const scalarField& volumes=subTimeMeshes[K].V();
            forAll(source,celli)
            {
                const label flat=phaseOffsets[K]+celli;
                phaseWeights[flat]=volumes[celli]/phaseVolume;
                phaseRhs[flat]=phaseIntegratedRhs[K][celli]/volumes[celli];
                const scalar difference=phaseIntegratedRhs[K][celli]-phaseExpectedRhs[K][celli];
                const scalar scale=phaseRhsScale[K][celli]+mag(source[celli])
                    +mag(phaseMatrixAction[K][celli]);
                assemblyFinite=assemblyFinite && std::isfinite(difference)
                    && std::isfinite(scale) && std::isfinite(phaseRhs[flat]);
                assemblyBackwardError=max(assemblyBackwardError,
                    mag(difference)/max(scale,scalar(VSMALL)));
                assemblyDifferenceSquare+=phaseWeights[flat]*sqr(difference/volumes[celli]);
                assemblyReferenceSquare+=phaseWeights[flat]
                    *sqr(phaseExpectedRhs[K][celli]/volumes[celli]);
            }
        }
        reduce(assemblyFinite,andOp<bool>());
        reduce(assemblyBackwardError,maxOp<scalar>());
        reduce(assemblyDifferenceSquare,sumOp<scalar>());
        reduce(assemblyReferenceSquare,sumOp<scalar>());
        const scalar assemblyRelative=Foam::sqrt(assemblyDifferenceSquare)
            /max(Foam::sqrt(assemblyReferenceSquare),scalar(VSMALL));
        Info<< "IHBPhaseTurbulenceAssembly field " << phaseFieldName
            << " nativeVsBlockBackwardError " << assemblyBackwardError
            << " relativeRmsDifference " << assemblyRelative
            << " finite " << assemblyFinite << endl;
        if (!assemblyFinite || !std::isfinite(assemblyRelative)
         || assemblyBackwardError>1e-10)
            FatalErrorInFunction << "phaseImplicit native/block RHS check failed for "
                << phaseFieldName << exit(FatalError);

        scalarField phaseIncrement(phaseFlatSize,Zero);
        bool phaseLinearConverged=false;
        scalar phaseLinearRelative=GREAT;
        #include "phaseTurbulenceFGMRES.H"

        // Global positivity protection acts on field increments, not equation
        // diagonals. It never clips a negative candidate into an accepted field.
        scalar positivityLimit=GREAT;
        forAll(subTimeMeshes,K)
        {
            const scalar minimum=hbTurbulence[K].phaseImplicitMinimum(phaseSolveEpsilon);
            forAll(phaseBase[K],celli)
            {
                const scalar dq=phaseIncrement[phaseOffsets[K]+celli];
                if (dq<0)
                    positivityLimit=min(positivityLimit,
                        (phaseBase[K][celli]-minimum)/(-dq));
            }
        }
        reduce(positivityLimit,minOp<scalar>());
        scalar damping=phaseTurbulenceDamping;
        if (positivityLimit<=damping) damping=0.99*positivityLimit;
        if (!std::isfinite(damping) || damping<phaseTurbulenceMinDamping)
        {
            if (phaseTurbulenceLineSearch)
            {
                phaseProposalValid=false;
                Info<< "IHBPhaseTurbulenceProposalRejected field " << phaseFieldName
                    << " reason positivityDampingBelowMinimum damping " << damping << endl;
                break;
            }
            FatalErrorInFunction
                << "phaseImplicit " << phaseFieldName
                << " positivity requires damping " << damping
                << " below phaseTurbulenceMinDamping " << phaseTurbulenceMinDamping
                << "; no clipped candidate has been committed" << exit(FatalError);
        }
        bool candidateValid=true;
        scalar maximumIncrement=0;
        forAll(subTimeMeshes,K)
        {
            const scalar minimum=hbTurbulence[K].phaseImplicitMinimum(phaseSolveEpsilon);
            forAll(phaseBase[K],celli)
            {
                const scalar increment=damping*phaseIncrement[phaseOffsets[K]+celli];
                const scalar candidate=phaseBase[K][celli]+increment;
                candidateValid=candidateValid && std::isfinite(candidate) && candidate>=minimum;
                maximumIncrement=max(maximumIncrement,mag(increment));
            }
        }
        reduce(candidateValid,andOp<bool>());
        reduce(maximumIncrement,maxOp<scalar>());
        if (!candidateValid)
        {
            if (phaseTurbulenceLineSearch)
            {
                phaseProposalValid=false;
                Info<< "IHBPhaseTurbulenceProposalRejected field " << phaseFieldName
                    << " reason invalidFieldCandidate" << endl;
                break;
            }
            FatalErrorInFunction << "Invalid phaseImplicit " << phaseFieldName
                << " candidate; field update rejected" << exit(FatalError);
        }
        forAll(subTimeMeshes,K)
        {
            scalarField& values=hbTurbulence[K].phaseImplicitField(phaseSolveEpsilon).primitiveFieldRef();
            forAll(values,celli)
                values[celli]=phaseBase[K][celli]+damping*phaseIncrement[phaseOffsets[K]+celli];
        }
        forAll(subTimeMeshes,K)
        {
            const HBPhaseTime phaseClock(runTime,snapshots[K]+hbMotionDeltaT,
                hbMotionDeltaT,subTimeMeshes[K],meshPhiHB[K]);
            volScalarField& q=hbTurbulence[K].phaseImplicitField(phaseSolveEpsilon);
            q.correctBoundaryConditions();
            const scalar minimum=hbTurbulence[K].phaseImplicitMinimum(phaseSolveEpsilon);
            forAll(q,celli)
                candidateValid=candidateValid && std::isfinite(q[celli]) && q[celli]>=minimum;
        }
        reduce(candidateValid,andOp<bool>());
        if (!candidateValid)
        {
            if (phaseTurbulenceLineSearch)
            {
                phaseProposalValid=false;
                Info<< "IHBPhaseTurbulenceProposalRejected field " << phaseFieldName
                    << " reason invalidBoundaryCandidate" << endl;
                break;
            }
            FatalErrorInFunction << "phaseImplicit boundary update produced invalid "
                << phaseFieldName << exit(FatalError);
        }
        Info<< "IHBPhaseTurbulenceUpdate field " << phaseFieldName
            << " damping " << damping << " configuredDamping " << phaseTurbulenceDamping
            << " positivityReduced " << (damping<phaseTurbulenceDamping)
            << " maxIncrement " << maximumIncrement
            << " linearConverged " << phaseLinearConverged
            << " trueRelative " << phaseLinearRelative
            << " nativeEquationRelaxation 0 timeRate 0"
            << " sharedPseudoTime " << phaseTurbulencePseudoTime
            << " proposalOnly " << phaseTurbulenceLineSearch
            << " linearResidualAppliesToUndampedDirection 1" << endl;
    }
    if (phaseProposalValid)
        forAll(subTimeMeshes,K)
        {
            const HBPhaseTime phaseClock(runTime,snapshots[K]+hbMotionDeltaT,
                hbMotionDeltaT,subTimeMeshes[K],meshPhiHB[K]);
            const bool valid=hbTurbulence[K].finishPhaseImplicit(!phaseTurbulenceLineSearch);
            phaseProposalValid=phaseProposalValid && valid;
        }
    phaseAcceptTurbulenceProposal();
}
