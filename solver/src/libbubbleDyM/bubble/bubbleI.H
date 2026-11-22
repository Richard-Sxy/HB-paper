/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::bubble::trackingData::trackingData
(
    const bubbleCloud& spc,
    const interpolationCellPoint<scalar>& rhoInterp,
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCellPoint<scalar>& nuInterp,
    const interpolationCellPoint<scalar>& pInterp,
    const interpolationCellPoint<scalar>& alphaCInterp,
    const interpolationCellPoint<vector>& curlUcInterp,
    const interpolationCellPoint<vector>& DUcDtInterp,
    const interpolationCellPoint<vector>& GradPcInterp,
    const vector& g
)
:
    particle::trackingData(spc),
    rhoInterp_(rhoInterp),
    UInterp_(UInterp),
    nuInterp_(nuInterp),
    pInterp_(pInterp),
    alphaCInterp_(alphaCInterp),
    curlUcInterp_(curlUcInterp),
    DUcDtInterp_(DUcDtInterp),
    GradPcInterp_(GradPcInterp),
    g_(g)
{}

inline Foam::bubble::RPODE::RPODE
(
    scalar p,
    scalar p0,
    scalar rho,
    scalar mu,
    scalar pSat,
    scalar sigma,
    scalar R0
)
:
   p_(p),
   p0_(p0),
   rho_(rho),
   mu_(mu_),
   pv_(pSat),
   sigma_(sigma),
   R0_(R0)
{}


inline Foam::bubble::bubble
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    particle(mesh, position, celli),
    d_(0.0),
    U_(Zero),
    U0_(Zero),
    R_(d_/2.0),
    R0_(R_),
    dRt_(0.0),
    ddRt_(0.0)

{}


inline Foam::bubble::bubble
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const scalar d,
    const vector& U
)
:
    particle(mesh, position, celli, tetFacei, tetPti),
    d_(d),
    U_(U),
    R_(d/2.0),
    R0_(R_),
    U0_(U),
    dRt_(0.0),
    ddRt_(0.0)
{}

inline Foam::bubble::bubble
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const scalar d,
    const scalar R0,
    const vector& U
)
:
    particle(mesh, position, celli),
    d_(d),
    U_(U),
    U0_(U),
    R_(d/2.0),
    R0_(R0),
    dRt_(0.0),
    ddRt_(0.0)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::bubble::trackingData::rhoInterp() const
{
    return rhoInterp_;
}


inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::bubble::trackingData::UInterp() const
{
    return UInterp_;
}


inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::bubble::trackingData::nuInterp() const
{
    return nuInterp_;
}

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::bubble::trackingData::pInterp() const
{
return pInterp_;
}

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::bubble::trackingData::alphaCInterp() const
{
	return alphaCInterp_;
}

inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::bubble::trackingData::curlUcInterp() const
{
    return curlUcInterp_;
}

inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::bubble::trackingData::DUcDtInterp() const
{
    return DUcDtInterp_;
}

inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::bubble::trackingData::GradPcInterp() const
{
    return GradPcInterp_;
}

inline const Foam::vector& Foam::bubble::trackingData::g() const
{
    return g_;
}

inline Foam::label Foam::bubble::RPODE::nEqns() const
{
     return 2;
}

inline void Foam::bubble::RPODE::derivatives
(
     const scalar x,
     const scalarField& y,
     scalarField& dydx
) const
{
    
     scalar kappa=4.0/3.0;
 
     if (y[0] >= R0_)
     {
          kappa=1;
     }
     scalar pB_ = (p0_ +2*sigma_/R0_-pv_)*Foam::pow(R0_/y[0],3*kappa)+pv_;
     dydx[0] = y[1];
     //dydx[1] = -17*pow(y[1],2)/(16*y[0]+VSMALL) + 2*(pB_- p_)/(rho_*y[0]+VSMALL) - 2*(2*sigma_ + 4*mu_*y[1])/(rho_*sqr(y[0]) + VSMALL);
     dydx[1] = -3*pow(y[1],2)/(2*y[0]+VSMALL) + (pB_- p_)/(rho_*y[0]+VSMALL) - (2*sigma_ + 4*mu_*y[1])/(rho_*sqr(y[0]) + VSMALL);
}
 
inline void Foam::bubble::RPODE::jacobian
(
     const scalar x,
     const scalarField& y,
     scalarField& dfdx,
     scalarSquareMatrix& dfdy
) const
{
     scalar kappa=4.0/3.0;

     if (y[0] >= R0_)
     {
         kappa=1;
     }

     scalar pB_ = (p0_+2*sigma_/R0_-pv_)*Foam::pow(R0_/y[0],(3*kappa))+pv_;
 
     dfdx[0] = 0.0;
     dfdx[1] = 0.0;

     dfdy(0, 0) = 0.0;
     dfdy(0, 1) = 1.0;
                                                                                                                                                           
     //dfdy(1, 0) = 17*pow(y[1],2)/(16*pow(y[0],2)+VSMALL) - 2*(pB_- p_)/(rho_*pow(y[0],2)+VSMALL) +2*(4*sigma_ + 8*mu_*y[1])/(rho_*pow(y[0],3) + VSMALL);
     //dfdy(1, 1) = -17*y[1]/(8*y[0]+VSMALL) - 8*mu_ /(rho_ * sqr(y[0])+VSMALL);
     dfdy(1, 0) = 3*pow(y[1],2)/(2*pow(y[0],2)+VSMALL) - (pB_- p_)/(rho_*pow(y[0],2)+VSMALL) +(4*sigma_ + 8*mu_*y[1])/(rho_*pow(y[0],3) + VSMALL);
     dfdy(1, 1) = -3*y[1]/(y[0]+VSMALL) - 4*mu_ /(rho_ * sqr(y[0])+VSMALL);
}


inline Foam::scalar& Foam::bubble::d()
{
    return d_;
}


inline Foam::vector& Foam::bubble::U()
{
    return U_;
}

inline Foam::vector& Foam::bubble::U0()
{
    return U0_;
}
inline Foam::scalar Foam::bubble::R()
{
    return R_;
}
inline Foam::scalar Foam::bubble::R0()
{
    return R0_;
}
inline Foam::scalar Foam::bubble::dRt()
{
    return dRt_;
}
inline Foam::scalar Foam::bubble::ddRt()
{
    return ddRt_;
}
// ************************************************************************* //
