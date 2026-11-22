/*---------------------------------------------------------------------------*\
    HB extension of OpenFOAM v2112 RASModels::kEpsilon.
    Native model equations: Copyright (C) 2011-2017 OpenFOAM Foundation,
    Copyright (C) 2019-2020 OpenCFD Ltd. GPL v3 or later, as ICSFoam.
\*---------------------------------------------------------------------------*/
#ifndef incompressibleHBKEpsilon_H
#define incompressibleHBKEpsilon_H
#include "kEpsilon.H"
#include "fvOptions.H"
#include "bound.H"
#include "IOHBZoneList.H"
#include "epsilonWallFunctionFvPatchScalarField.H"
#include "cyclicAMIFvPatch.H"
#include "processorFvPatch.H"

namespace Foam
{
// Preserve native kEpsilon coefficients, fields, wall functions and model name.
// Phase coupling uses synchronised three-stage pseudo-time sweeps (SSP RK3).
class HBKEpsilon : public RASModels::kEpsilon<incompressible::turbulenceModel>
{
    typedef incompressible::turbulenceModel Basic;
    typedef RASModels::kEpsilon<Basic> Base;
    autoPtr<volScalarField::Internal> GbyNu_, G_;
    scalarField kTimeSource_, epsTimeSource_, timeDiag_, timeRate_;
    scalarField kFrozen_, epsFrozen_;

    tmp<fvScalarMatrix> equation(const bool eps)
    {
        fv::options& options = fv::options::New(this->mesh_);
        const volScalarField::Internal divU
        (
            fvc::div(fvc::absolute(this->phi(), this->U_))().v()
        );
        // Unit alpha/rho: incompressible single-phase model only.
        // Physical d/dt is supplied below by HB, never by pseudo runTime.
        if (eps)
            return
            (
                fvm::div(this->alphaRhoPhi_, epsilon_)
              - fvm::laplacian(DepsilonEff(), epsilon_)
             == C1_*GbyNu_()*Cmu_*k_()
              - fvm::SuSp(((2.0/3.0)*C1_-C3_)*divU, epsilon_)
              - fvm::Sp(C2_*epsilon_()/k_(), epsilon_)
              + epsilonSource() + options(this->alpha_, this->rho_, epsilon_)
            );
        return
        (
            fvm::div(this->alphaRhoPhi_, k_)
          - fvm::laplacian(DkEff(), k_)
         == G_()
          - fvm::SuSp((2.0/3.0)*divU, k_)
          - fvm::Sp(epsilon_()/k_(), k_)
          + kSource() + options(this->alpha_, this->rho_, k_)
        );
    }

    void addTime(fvScalarMatrix& eqn, const bool eps, const bool stabilise)
    {
        const scalarField& frozen = eps ? epsFrozen_ : kFrozen_;
        eqn.source() += (eps ? epsTimeSource_ : kTimeSource_);
        // Move D_jj q_j from frozen RHS to the implicit diagonal.
        eqn.diag() += timeDiag_;
        eqn.source() += timeDiag_*frozen;
        if (stabilise)
        {
            // Cancels exactly at the fixed point; excluded from residuals.
            eqn.diag() += timeRate_;
            eqn.source() += timeRate_*frozen;
        }
    }

    void constrain(fvScalarMatrix& eqn, const bool eps)
    {
        fv::options::New(this->mesh_).constrain(eqn);
        // Must follow the HB addition: native prescribed wall rows replace
        // the COMPLETE equation, including its temporal source/diagonal.
        if (eps) eqn.boundaryManipulate(epsilon_.boundaryFieldRef());
    }


    // Shared assembly for ordinary diagnostics and candidate acceptance.
    // freezeTime must have been refreshed for ALL phases before this call.
    scalarField physicalResidualField(const bool eps)
    {
        tmp<fvScalarMatrix> tq=equation(eps);
        addTime(tq.ref(),eps,false);
        constrain(tq.ref(),eps);
        const scalarField& q=eps ? epsilon_.primitiveField() : k_.primitiveField();
        const volScalarField H(tq().H());
        const volScalarField A(tq().A());
        return scalarField(H.primitiveField()-A.primitiveField()*q);
    }

    static void restorePhaseField(volScalarField& live,const volScalarField& saved)
    {
        live.primitiveFieldRef()=saved.primitiveField();
        auto& patches=live.boundaryFieldRef();
        forAll(patches,patchi)
            patches.set(patchi,saved.boundaryField()[patchi].clone(live.internalField()).ptr());
        // Native clone preserves BC coefficients/values while resetting update
        // and matrix-manipulation flags; epsilon wall caches are rebuilt on use.
    }

    // Diagnostics only: reuse the already assembled, unrelaxed residual.
    // No boundary updates, matrix construction or field writes occur here.
    void auditEpsilonResidual(const scalarField& residual, const label phase) const
    {
        const scalarField& volumes = this->mesh_.V();
        // Initialise any geometry cache on every rank before selecting a winner.
        const vectorField& centres = this->mesh_.C().primitiveField();
        scalar localPeak = -1;
        label localCell = -1;
        scalarList data(15, Zero);
        forAll(residual, celli)
        {
            const scalar magnitude = mag(residual[celli]);
            if (magnitude > localPeak)
            {
                localPeak = magnitude;
                localCell = celli;
            }
            data[13] += volumes[celli];
            data[14] += volumes[celli]*sqr(residual[celli]);
        }
        scalar globalPeak = localPeak;
        reduce(globalPeak, maxOp<scalar>());
        label winner = localCell >= 0 && localPeak == globalPeak
            ? Pstream::myProcNo() : Pstream::nProcs();
        reduce(winner, minOp<label>());

        // A tie selects the smallest rank, then the smallest local cell.
        // All ranks still participate in the same collective operations.
        if (Pstream::myProcNo() == winner)
        {
            const vector& centre = centres[localCell];
            data[0] = localCell;
            data[1] = centre.x();
            data[2] = centre.y();
            data[3] = centre.z();
            data[4] = volumes[localCell];
            data[5] = residual[localCell];
            data[6] = epsilon_.primitiveField()[localCell];
            data[7] = k_.primitiveField()[localCell];
            data[8] = this->nut_.primitiveField()[localCell];
            data[9] = epsTimeSource_[localCell]/volumes[localCell];
            forAll(this->mesh_.boundary(), patchi)
            {
                const fvPatch& patch = this->mesh_.boundary()[patchi];
                for (const label celli : patch.faceCells())
                {
                    if (celli != localCell) continue;
                    if (isA<epsilonWallFunctionFvPatchScalarField>
                        (epsilon_.boundaryField()[patchi])) data[10] = 1;
                    if (isA<cyclicAMIFvPatch>(patch)) data[11] = 1;
                    if (isA<processorFvPatch>(patch)) data[12] = 1;
                    break;
                }
            }
        }
        Pstream::listCombineGather(data, plusEqOp<scalar>());
        Pstream::listCombineScatter(data);
        if (globalPeak < 0 || data[13] <= 0)
        {
            Info<< "IHBEpsilonResidualAudit phase " << phase
                << " emptyMesh 1" << endl;
            return;
        }
        // epsilon wall-function rows replace the complete HB equation.
        const scalar hbResidual = data[10] > 0 ? 0 : data[9];
        Info<< "IHBEpsilonResidualAudit phase " << phase
            << " rank " << winner << " cell " << label(data[0])
            << " C " << vector(data[1], data[2], data[3])
            << " V " << data[4] << " maxAbs " << globalPeak
            << " signedResidual " << data[5]
            << " epsilon " << data[6] << " k " << data[7]
            << " nut " << data[8]
            << " hbResidual " << hbResidual
            << " spatialResidual " << data[5]-hbResidual
            << " unconstrainedHbResidual " << data[9]
            << " wallConstraint " << label(data[10])
            << " adjacentAMI " << label(data[11])
            << " adjacentProcessor " << label(data[12])
            << " volumeRms " << sqrt(data[14]/data[13]) << endl;
    }

public:
    HBKEpsilon(const volVectorField& U, const surfaceScalarField& phi,
               const transportModel& transport)
    : Base(geometricOneField(), geometricOneField(), U, phi, phi, transport),
      kTimeSource_(U.size(), Zero), epsTimeSource_(U.size(), Zero),
      timeDiag_(U.size(), Zero), timeRate_(U.size(), Zero),
      kFrozen_(U.size(), Zero), epsFrozen_(U.size(), Zero)
    {}

    bool active() const { return this->turbulence_; }

    void prepare()
    {
        eddyViscosity<RASModel<Basic>>::correct();
        k_.correctBoundaryConditions();
        G_.clear(); GbyNu_.clear();
        const tmp<volTensorField> gradU = fvc::grad(this->U_);
        GbyNu_.reset(new volScalarField::Internal
        (
            this->type()+":GbyNu",
            gradU().v() && dev(twoSymm(gradU().v()))
        ));
        G_.reset(new volScalarField::Internal
        (
            this->GName(), this->nut_()*GbyNu_()
        ));
        // G must be registered while epsilon wall functions update it.
        epsilon_.correctBoundaryConditions();
        epsilon_.boundaryFieldRef().updateCoeffs();
    }

    // Call for ALL phases before solving ANY phase (Jacobi phase snapshot).
    void freezeTime(const HBZoneList& HB, const UPtrList<HBKEpsilon>& models,
                    const label row, const scalar damping)
    {
        kTimeSource_=Zero; epsTimeSource_=Zero;
        timeDiag_=Zero; timeRate_=Zero;
        kFrozen_=k_.primitiveField(); epsFrozen_=epsilon_.primitiveField();
        forAll(models, col)
        {
            blockFvMatrix<scalar,scalar> block(this->mesh_);
            block.diag()=Zero;
            HB.addBlock(block,row,col); // V*D, including per-zone selection
            const scalarField& vd=block.diag();
            kTimeSource_ -= vd*models[col].k_.primitiveField();
            epsTimeSource_ -= vd*models[col].epsilon_.primitiveField();
            timeRate_ += damping*mag(vd);
            if (col==row) timeDiag_=vd;
        }
    }

    virtual void correct() override
    {
        if (!active()) return;
        fv::options& options=fv::options::New(this->mesh_);
        for (const bool eps : {true,false})
        {
            volScalarField& q=eps ? epsilon_ : k_;
            tmp<fvScalarMatrix> tq=equation(eps);
            fvScalarMatrix& eqn=tq.ref();
            addTime(eqn,eps,true);
            eqn.relax();
            constrain(eqn,eps);
            eqn.solve();
            options.correct(q);
            bound(q,eps ? this->epsilonMin_ : this->kMin_);
        }
        correctNut();
    }

    // Full numeric transaction state, including boundary coefficients and the
    // cached HB/G data. Saved fields are unregistered and are never evaluated.
    struct PhaseImplicitState
    {
        volScalarField k,epsilon,nut;
        scalarField gByNu,g,kTimeSource,epsTimeSource,timeDiag,timeRate,kFrozen,epsFrozen;
        explicit PhaseImplicitState(const HBKEpsilon& model)
        :
            k(IOobject("phaseSavedK",model.mesh_.time().timeName(),model.mesh_,
                IOobject::NO_READ,IOobject::NO_WRITE,false),model.k_),
            epsilon(IOobject("phaseSavedEpsilon",model.mesh_.time().timeName(),model.mesh_,
                IOobject::NO_READ,IOobject::NO_WRITE,false),model.epsilon_),
            nut(IOobject("phaseSavedNut",model.mesh_.time().timeName(),model.mesh_,
                IOobject::NO_READ,IOobject::NO_WRITE,false),model.nut_),
            gByNu(model.GbyNu_().field()),g(model.G_().field()),
            kTimeSource(model.kTimeSource_),epsTimeSource(model.epsTimeSource_),
            timeDiag(model.timeDiag_),timeRate(model.timeRate_),
            kFrozen(model.kFrozen_),epsFrozen(model.epsFrozen_)
        {}
    };

    void restorePhaseImplicitState(const PhaseImplicitState& saved)
    {
        restorePhaseField(k_,saved.k);
        restorePhaseField(epsilon_,saved.epsilon);
        restorePhaseField(this->nut_,saved.nut);
        GbyNu_().field()=saved.gByNu;
        G_().field()=saved.g;
        kTimeSource_=saved.kTimeSource; epsTimeSource_=saved.epsTimeSource;
        timeDiag_=saved.timeDiag; timeRate_=saved.timeRate;
        kFrozen_=saved.kFrozen; epsFrozen_=saved.epsFrozen;
    }

    // Local difference only; caller reduces once across phases/ranks. Includes
    // caches so lambda=0 cannot silently retain a candidate's frozen HB source.
    scalar phaseImplicitStateDifference
    (
        const PhaseImplicitState& saved,const bool includeCaches=true
    ) const
    {
        scalar maximum=0;
        const auto compare = [&](const scalarField& a,const scalarField& b)
        {
            forAll(a,i)
            {
                const scalar difference=mag(a[i]-b[i])/max(scalar(1),mag(b[i]));
                if (!std::isfinite(difference)) maximum=GREAT;
                else maximum=max(maximum,difference);
            }
        };
        const auto compareField = [&](const volScalarField& a,const volScalarField& b)
        {
            compare(a.primitiveField(),b.primitiveField());
            forAll(a.boundaryField(),patchi)
                compare(a.boundaryField()[patchi],b.boundaryField()[patchi]);
        };
        compareField(k_,saved.k); compareField(epsilon_,saved.epsilon);
        compareField(this->nut_,saved.nut);
        if (!includeCaches) return maximum;
        compare(GbyNu_().field(),saved.gByNu); compare(G_().field(),saved.g);
        compare(kTimeSource_,saved.kTimeSource); compare(epsTimeSource_,saved.epsTimeSource);
        compare(timeDiag_,saved.timeDiag); compare(timeRate_,saved.timeRate);
        compare(kFrozen_,saved.kFrozen); compare(epsFrozen_,saved.epsFrozen);
        return maximum;
    }

    const volScalarField& phaseImplicitNut() const { return this->nut_; }

    // A trial is rebuilt from the same saved fields before entering here.
    // Unlike finishPhaseImplicit(), invalid trials return false for rollback.
    bool preparePhaseImplicitResidual()
    {
        bool valid=true;
        for (const bool eps : {true,false})
        {
            const scalarField& q=eps ? epsilon_.primitiveField() : k_.primitiveField();
            const scalar minimum=eps ? this->epsilonMin_.value() : this->kMin_.value();
            for (const scalar value : q)
                valid=valid && std::isfinite(value) && value>=minimum;
        }
        reduce(valid,andOp<bool>());
        if (!valid) return false;
        k_.correctBoundaryConditions();
        epsilon_.correctBoundaryConditions();
        correctNut();
        for (const scalar value : this->nut_.primitiveField())
            valid=valid && std::isfinite(value) && value>=0;
        for (const fvPatchScalarField& patch : this->nut_.boundaryField())
            for (const scalar value : patch)
                valid=valid && std::isfinite(value) && value>=0;
        reduce(valid,andOp<bool>());
        if (!valid) return false;
        // Rebuild G from this NEW nut, then reapply the native wall production.
        prepare();
        for (const bool eps : {true,false})
        {
            const scalarField& q=eps ? epsilon_.primitiveField() : k_.primitiveField();
            const scalar minimum=eps ? this->epsilonMin_.value() : this->kMin_.value();
            for (const scalar value : q)
                valid=valid && std::isfinite(value) && value>=minimum;
        }
        for (const scalar value : G_().field()) valid=valid && std::isfinite(value);
        for (const scalar value : GbyNu_().field()) valid=valid && std::isfinite(value);
        reduce(valid,andOp<bool>());
        return valid;
    }

    // Add LOCAL integral statistics; the caller performs all MPI reductions.
    bool phaseImplicitResidualStatistics
    (
        scalar& kSquare,scalar& epsSquare,scalar& volume,
        scalarField* kResidual=nullptr,scalarField* epsilonResidual=nullptr
    )
    {
        bool finite=true;
        for (const scalar cellVolume : this->mesh_.V()) volume+=cellVolume;
        for (const bool eps : {true,false})
        {
            const scalarField residual(physicalResidualField(eps));
            scalarField* saved=eps ? epsilonResidual : kResidual;
            if (saved) *saved=residual;
            const scalarField& q=eps ? epsilon_.primitiveField() : k_.primitiveField();
            const scalar minimum=eps ? this->epsilonMin_.value() : this->kMin_.value();
            scalar& square=eps ? epsSquare : kSquare;
            forAll(residual,celli)
            {
                finite=finite && std::isfinite(residual[celli])
                    && std::isfinite(q[celli]) && q[celli]>=minimum;
                square+=this->mesh_.V()[celli]*sqr(residual[celli]);
            }
        }
        return finite && std::isfinite(kSquare) && std::isfinite(epsSquare)
            && std::isfinite(volume);
    }

    // Narrow API for the optional all-phase scalar solve. This deliberately
    // bypasses native equation relaxation and the SSPRK pseudo-time diagonal.
    volScalarField& phaseImplicitField(const bool eps)
    {
        return eps ? epsilon_ : k_;
    }

    scalar phaseImplicitMinimum(const bool eps) const
    {
        return eps ? this->epsilonMin_.value() : this->kMin_.value();
    }

    tmp<fvScalarMatrix> phaseImplicitEquation(const bool eps)
    {
        return equation(eps);
    }

    void phaseImplicitConstrain(fvScalarMatrix& eqn, const bool eps)
    {
        constrain(eqn,eps);
    }

    // No generic inference of constrained cells from zero matrix offdiagonals:
    // only the native, unweighted epsilonWallFunction setValues is supported.
    boolList phaseImplicitConstraintMask(const bool eps) const
    {
        boolList mask(this->mesh_.nCells(),false);
        if (eps)
            forAll(epsilon_.boundaryField(),patchi)
                if (epsilon_.boundaryField()[patchi].type()=="epsilonWallFunction")
                    for (const label celli : this->mesh_.boundary()[patchi].faceCells())
                        mask[celli]=true;
        return mask;
    }

    void validatePhaseImplicit() const
    {
        if (!active())
            FatalErrorInFunction << "phaseImplicit requires active kEpsilon in every phase"
                << exit(FatalError);
        if (!static_cast<const fv::optionList&>
            (fv::options::New(this->mesh_)).empty())
            FatalErrorInFunction
                << "phaseImplicit currently requires an empty fvOptions list; "
                << "unknown cell constraints cannot be inferred safely"
                << exit(FatalError);
        for (const bool eps : {true,false})
        {
            const volScalarField& q=eps ? epsilon_ : k_;
            forAll(q.boundaryField(),patchi)
            {
                const word type=q.boundaryField()[patchi].type();
                const bool supported=type=="fixedValue" || type=="zeroGradient"
                    || type=="inletOutlet" || type=="empty"
                    || type=="cyclicAMI" || type=="processor"
                    || (eps && type=="epsilonWallFunction")
                    || (!eps && type=="kqRWallFunction");
                if (!supported)
                    FatalErrorInFunction
                        << "Unsupported phaseImplicit boundary " << q.name()
                        << '/' << q.boundaryField()[patchi].patch().name()
                        << " of type " << type << "; constraint handling must "
                        << "be reviewed before extending this experimental mode"
                        << exit(FatalError);
            }
        }
    }

    bool finishPhaseImplicit(const bool fatalOnInvalid=true)
    {
        // k was updated after the epsilon solve. Recompute prescribed epsilon
        // wall values from this new k BEFORE deriving nut, exactly as blend()
        // does. G_ remains registered, and the caller holds this phase's clock.
        epsilon_.correctBoundaryConditions();
        bool positive=true;
        for (const scalar value : epsilon_.primitiveField())
            positive=positive && std::isfinite(value)
                && value>=this->epsilonMin_.value();
        reduce(positive,andOp<bool>());
        if (!positive)
        {
            if (!fatalOnInvalid) return false;
            FatalErrorInFunction
                << "phaseImplicit wall refresh produced invalid epsilon"
                << exit(FatalError);
        }
        correctNut();
        bool valid=true;
        for (const scalar value : this->nut_.primitiveField())
            valid=valid && std::isfinite(value) && value>=0;
        for (const fvPatchScalarField& patch : this->nut_.boundaryField())
            for (const scalar value : patch)
                valid=valid && std::isfinite(value) && value>=0;
        reduce(valid,andOp<bool>());
        if (!valid)
        {
            if (!fatalOnInvalid) return false;
            FatalErrorInFunction
                << "phaseImplicit produced invalid nut after positive k/epsilon update"
                << exit(FatalError);
        }
        return true;
    }

    const scalarField& kValues() const { return k_.primitiveField(); }
    const scalarField& epsilonValues() const { return epsilon_.primitiveField(); }
    const scalarField& dissipationValues() const { return epsilonValues(); }

    void blend(const scalar weight, const scalarField& kStart,
               const scalarField& epsStart)
    {
        k_.primitiveFieldRef()=weight*kStart+(1-weight)*k_.primitiveField();
        epsilon_.primitiveFieldRef()=weight*epsStart+(1-weight)*epsilon_.primitiveField();
        k_.correctBoundaryConditions();
        epsilon_.correctBoundaryConditions();
        correctNut();
    }

    // Unrelaxed physical equations, with current all-phase HB sources.
    // H-A*q includes coupled neighbour values ONCE, just like momentum.
    // v2112 fvScalarMatrix::residual also adds coupled boundary sources after
    // its ldu interface action; it is not used for this physical residual.
    void residuals
    (
        scalar& kMax, scalar& epsMax,
        const bool audit = false, const label phase = -1
    )
    {
        for (const bool eps : {true,false})
        {
            const scalarField residual(physicalResidualField(eps));
            const scalarField& q=eps ? epsilon_.primitiveField() : k_.primitiveField();
            scalar& result=eps ? epsMax : kMax;
            forAll(residual,c)
            {
                const scalar error=mag(residual[c]);
                if (!std::isfinite(error) || !std::isfinite(q[c]) || q[c]<=0)
                    FatalErrorInFunction << "Invalid HB turbulence field/residual: "
                        << (eps ? "epsilon" : "k") << exit(FatalError);
                result=max(result,error);
            }
            if (eps && audit) auditEpsilonResidual(residual, phase);
        }
    }
};
}
#endif
