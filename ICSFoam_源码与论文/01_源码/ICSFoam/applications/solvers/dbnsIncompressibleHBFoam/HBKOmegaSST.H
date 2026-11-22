/*---------------------------------------------------------------------------*\
    Harmonic-balance extension of OpenCFD OpenFOAM v2112 kOmegaSST.
    Spatial equations derived from kOmegaSSTBase.C:
    Copyright (C) 2011-2016 OpenFOAM Foundation, 2017-2021 OpenCFD Ltd.
    GPL v3 or later, as ICSFoam.
\*---------------------------------------------------------------------------*/
#ifndef incompressibleHBKOmegaSST_H
#define incompressibleHBKOmegaSST_H

#include "kOmegaSST.H"
#include "fvOptions.H"
#include "bound.H"
#include "IOHBZoneList.H"

namespace Foam
{

// Single-phase, constant-density SST. Native constitutive functions and wall
// functions are inherited; physical d/dt is replaced by the HB operator.
class HBKOmegaSST
:
    public RASModels::kOmegaSST<incompressible::turbulenceModel>
{
    typedef incompressible::turbulenceModel Basic;
    typedef RASModels::kOmegaSST<Basic> Base;

    autoPtr<volTensorField> gradUHB_;
    autoPtr<volScalarField> S2HB_, CDHB_, F1HB_, F23HB_;
    autoPtr<volScalarField::Internal> GbyNuHB_, GHB_;
    scalarField kTimeSource_, omegaTimeSource_, timeDiag_, timeRate_;
    scalarField kFrozen_, omegaFrozen_;

    tmp<fvScalarMatrix> equation(const bool frequency)
    {
        fv::options& options=fv::options::New(this->mesh_);
        const volScalarField::Internal divU
        (
            fvc::div(fvc::absolute(this->phi(), this->U_))().v()
        );
        if (frequency)
        {
            const volScalarField::Internal gamma(this->gamma(F1HB_().v()));
            const volScalarField::Internal beta(this->beta(F1HB_().v()));
            const volScalarField::Internal limitedGbyNu
            (
                this->GbyNu(GbyNuHB_(), F23HB_().v(), S2HB_().v())
            );
            return
            (
                fvm::div(this->alphaRhoPhi_, this->omega_)
              - fvm::laplacian(this->DomegaEff(F1HB_()), this->omega_)
             == gamma*limitedGbyNu
              - fvm::SuSp((2.0/3.0)*gamma*divU, this->omega_)
              - fvm::Sp(beta*this->omega_(), this->omega_)
              - fvm::SuSp
                (
                    (F1HB_().v()-scalar(1))*CDHB_().v()/this->omega_(),
                    this->omega_
                )
              + beta*sqr(this->omegaInf_)
              + this->Qsas(S2HB_().v(), gamma, beta)
              + this->omegaSource()
              + options(this->alpha_, this->rho_, this->omega_)
            );
        }
        return
        (
            fvm::div(this->alphaRhoPhi_, this->k_)
          - fvm::laplacian(this->DkEff(F1HB_()), this->k_)
         == this->Pk(GHB_())
          - fvm::SuSp((2.0/3.0)*divU, this->k_)
          - fvm::Sp(this->epsilonByk(F1HB_(), gradUHB_()), this->k_)
          + this->betaStar_*this->omegaInf_*this->kInf_
          + this->kSource()
          + options(this->alpha_, this->rho_, this->k_)
        );
    }

    void addTime(fvScalarMatrix& eqn, const bool frequency, const bool stabilise)
    {
        const scalarField& frozen=frequency ? omegaFrozen_ : kFrozen_;
        eqn.source() += frequency ? omegaTimeSource_ : kTimeSource_;
        // freezeTime stores -sum_l V*D_jl*q_l. Treat the self term implicitly.
        eqn.diag() += timeDiag_;
        eqn.source() += timeDiag_*frozen;
        if (stabilise)
        {
            // Pseudo-time regularisation cancels at the fixed point.
            eqn.diag() += timeRate_;
            eqn.source() += timeRate_*frozen;
        }
    }

    void constrain(fvScalarMatrix& eqn, const bool frequency)
    {
        fv::options::New(this->mesh_).constrain(eqn);
        // Native omega wall rows replace the COMPLETE equation, including HB.
        if (frequency) eqn.boundaryManipulate(this->omega_.boundaryFieldRef());
    }

public:
    HBKOmegaSST
    (
        const volVectorField& U, const surfaceScalarField& phi,
        const transportModel& transport
    )
    :
        Base(geometricOneField(), geometricOneField(), U, phi, phi, transport),
        kTimeSource_(U.size(), Zero), omegaTimeSource_(U.size(), Zero),
        timeDiag_(U.size(), Zero), timeRate_(U.size(), Zero),
        kFrozen_(U.size(), Zero), omegaFrozen_(U.size(), Zero)
    {}

    bool active() const { return this->turbulence_; }
    const scalarField& kValues() const { return this->k_.primitiveField(); }
    const scalarField& dissipationValues() const
    {
        return this->omega_.primitiveField();
    }

    // Same coefficient/production ordering as native SST. Keep G registered
    // while omegaWallFunction updates production and near-wall omega.
    void prepare()
    {
        eddyViscosity<RASModel<Basic>>::correct();
        this->k_.correctBoundaryConditions();
        GHB_.clear(); GbyNuHB_.clear();
        F1HB_.clear(); F23HB_.clear(); CDHB_.clear();
        S2HB_.clear(); gradUHB_.clear();
        gradUHB_.reset(new volTensorField("hbSST:gradU", fvc::grad(this->U_)));
        S2HB_.reset(new volScalarField
        (
            "hbSST:S2", 2*magSqr(symm(gradUHB_()))
        ));
        GbyNuHB_.reset(new volScalarField::Internal
        (
            this->type()+":GbyNu",
            gradUHB_().v() && dev(twoSymm(gradUHB_().v()))
        ));
        GHB_.reset(new volScalarField::Internal
        (
            this->GName(), this->nut_()*GbyNuHB_()
        ));
        this->omega_.correctBoundaryConditions();
        this->omega_.boundaryFieldRef().updateCoeffs();
        CDHB_.reset(new volScalarField
        (
            "hbSST:CDkOmega",
            (2*this->alphaOmega2_)
           *(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
        ));
        F1HB_.reset(new volScalarField("hbSST:F1", this->F1(CDHB_())));
        F23HB_.reset(new volScalarField("hbSST:F23", this->F23()));
    }

    // Snapshot ALL phases before solving any: no Gauss-Seidel time-source lag.
    void freezeTime
    (
        const HBZoneList& HB, const UPtrList<HBKOmegaSST>& models,
        const label row, const scalar damping
    )
    {
        kTimeSource_=Zero; omegaTimeSource_=Zero;
        timeDiag_=Zero; timeRate_=Zero;
        kFrozen_=kValues(); omegaFrozen_=dissipationValues();
        forAll(models,col)
        {
            blockFvMatrix<scalar,scalar> block(this->mesh_);
            block.diag()=Zero;
            HB.addBlock(block,row,col);
            const scalarField& vd=block.diag();
            kTimeSource_ -= vd*models[col].kValues();
            omegaTimeSource_ -= vd*models[col].dissipationValues();
            timeRate_ += damping*mag(vd);
            if (col==row) timeDiag_=vd;
        }
    }

    virtual void correct() override
    {
        if (!active()) return;
        fv::options& options=fv::options::New(this->mesh_);
        // Native SST solves omega before k, then refreshes the viscosity.
        for (const bool frequency : {true,false})
        {
            volScalarField& q=frequency ? this->omega_ : this->k_;
            tmp<fvScalarMatrix> tq=equation(frequency);
            addTime(tq.ref(),frequency,true);
            tq.ref().relax();
            constrain(tq.ref(),frequency);
            tq.ref().solve();
            options.correct(q);
            bound(q,frequency ? this->omegaMin_ : this->kMin_);
        }
        this->correctNut(S2HB_());
    }

    void blend
    (
        const scalar weight, const scalarField& kStart,
        const scalarField& omegaStart
    )
    {
        this->k_.primitiveFieldRef()=weight*kStart+(1-weight)*kValues();
        this->omega_.primitiveFieldRef()=
            weight*omegaStart+(1-weight)*dissipationValues();
        this->k_.correctBoundaryConditions();
        this->omega_.correctBoundaryConditions();
        this->correctNut(S2HB_());
    }

    // Unrelaxed physical residuals per unit volume: no pseudo-time damping.
    void residuals
    (
        scalar& kMax, scalar& omegaMax,
        const bool audit=false, const label phase=-1
    )
    {
        for (const bool frequency : {true,false})
        {
            tmp<fvScalarMatrix> tq=equation(frequency);
            addTime(tq.ref(),frequency,false);
            constrain(tq.ref(),frequency);
            const scalarField& q=frequency ? dissipationValues() : kValues();
            const volScalarField H(tq().H()), A(tq().A());
            const scalarField residual(H.primitiveField()-A.primitiveField()*q);
            scalar& peak=frequency ? omegaMax : kMax;
            forAll(residual,celli)
            {
                if (!std::isfinite(residual[celli]) || !std::isfinite(q[celli])
                 || q[celli]<=0)
                    FatalErrorInFunction << "Invalid SST-HB field/residual: "
                        << (frequency ? "omega" : "k") << exit(FatalError);
                peak=max(peak,mag(residual[celli]));
            }
        }
        if (audit)
            Info<< "IHBSSTResidualAudit phase " << phase
                << " kMin " << gMin(kValues())
                << " omegaMin " << gMin(dissipationValues())
                << " nutMin " << gMin(this->nut_.primitiveField()) << endl;
    }
};

}
#endif
