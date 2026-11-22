/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "bubbleCloud.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<bubble>, 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::bubble::move
(
    bubbleCloud& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
    {
        if (debug)
        {
            Pout<< "Time = " << mesh().time().timeName()
                << " trackTime = " << trackTime
                << " steptFraction() = " << stepFraction() << endl;
        }
        Pout<<"";
        const scalar sfrac = stepFraction();

        const scalar f = 1 - stepFraction();

        trackToAndHitFace(f*trackTime*U_, f, cloud, td);

        const scalar dt = (stepFraction() - sfrac)*trackTime;
        if(dt == 0)
        {
           return td.keepParticle;
        }
        const tetIndices tetIs = this->currentTetIndices();
        scalar rhoc = td.rhoInterp().interpolate(this->coordinates(), tetIs);
        //scalar rhoc = 998.85;
        vector Uc = td.UInterp().interpolate(this->coordinates(), tetIs);
        scalar nuc = td.nuInterp().interpolate(this->coordinates(), tetIs);
       // scalar nuc = 0.000001;
        scalar pc = td.pInterp().interpolate(this->coordinates(), tetIs);
	vector curlUc = td.curlUcInterp().interpolate(this->coordinates(), tetIs);

	vector DUcDt = td.DUcDtInterp().interpolate(this->coordinates(), tetIs);

	vector GradPc = td.GradPcInterp().interpolate(this->coordinates(), tetIs);
        
        vector posb=position();
        scalar pc_avg; 

        label numcell=0;
        for( int i = 0; i<6;i++)
        {
               int sign=pow(-1,double(i));
               int coord=i/2;
               vector pos2R = this->position();
               pos2R.component(coord) += sign * 2*R_ ;
               label cellI = -1;
               label tetFaceI = -1;
               label tetPtI = -1;

               mesh().findCellFacePt(pos2R, cellI, tetFaceI, tetPtI);
               if(cellI > -1)
               {
                  numcell += 1; 
                  pc_avg += td.pInterp().interpolate(pos2R, cellI);

               }

        }
        if( numcell > 0)
        {
             pc_avg = pc_avg/numcell;     
        }
        if(pc_avg == 0)
        {
            pc_avg = pc;
        }
        // bubble motion 
	scalar rhop = cloud.rhop();
        scalar p0 = cloud.p0();
        scalar pv = cloud.pSat();
        scalar sigma = cloud.sigma();
  	 
	    scalar magUr = mag(Uc - U_);

            scalar Re = magUr*d_/(nuc+VSMALL);
	    scalar mass = rhop*pi*pow(R_,3)*4.0/3.0;

	    scalar muc = nuc*rhoc;
        
            U0_ = U_;
            // Drag Force
	    scalar Dc = 0.75*muc*CdRe(Re)/(rhop*sqr(d_));
	    //Lift Force --?
            vector Lc = 6.46*rhoc*Cl(Re,rhoc,muc,curlUc)/rhop*(Uc-U_)^curlUc;
            scalar Vc = 3.0/(2.0*R_*rhop+VSMALL)*rhoc*dRt_;
        scalar dtX = dt / (1.0 + rhoc/(2.0 * rhop)+VSMALL);
        Info<<"";
        U_ = (U_ + dtX*(Dc*Uc + Lc  + rhoc/(2*rhop+VSMALL)*DUcDt + (1.0 - rhoc/(rhop+VSMALL))*td.g()+Vc*Uc  - GradPc/rhop))/(1.0 + dtX *(Vc + Dc) +VSMALL);
        if(mag(U_) < 1e-10)
        {
             U_ = vector::zero;
        } 
        
        bool solveRP = cloud.solveRP();
        // solver Rayleigh - Plesset Equation
        if(td.keepParticle == true && solveRP)
        {
         	#include "RP_Solver.H"
        }
    }
    return td.keepParticle;
}

// Drag Coefficient --(Amsden, O’Rourke & Butler 1989)
Foam::scalar Foam::bubble::CdRe
(
     scalar Re
)
{
     if(Re > 1000)
     {
	 return 0.424*Re;
     }
     else
     {
	 return 24.0*(1.0 + 1.0/6.0*pow(Re, 2.0/3.0));
     }
}

// Lift Coefficient  --(Mei,1992)
Foam::scalar Foam::bubble::Cl
(
    scalar Re,
    scalar rhoc,
    scalar muc,
    vector& curlUc
)
{
    scalar Rew = rhoc*mag(curlUc)*sqr(d_)/(muc+ VSMALL);
    scalar beta = 0.5*(Rew/(Re+ VSMALL));
    scalar alpha = 0.3314*sqrt(beta);
    scalar f = (1.0 - alpha)*exp(-0.1*Re)+alpha;
    if(Re > 40)
    {
        return 0.1112/(2*pi);
    }
    else
    {
	return 3.0/(2*pi*sqrt(Rew))*f; 
    }
}



bool Foam::bubble::hitPatch(bubbleCloud&, trackingData&)
{
    return false;
}


void Foam::bubble::hitProcessorPatch
(
    bubbleCloud&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::bubble::hitWallPatch(bubbleCloud& cloud, trackingData&)
{
    const vector nw = normal();

    scalar Un = U_ & nw;
    vector Ut = U_ - Un*nw;

    if (Un > 0)
    {
        U_ -= (1.0 + cloud.e())*Un*nw;
    }

    U_ -= cloud.mu()*Ut;
}


void Foam::bubble::transformProperties(const tensor& T)
{
    particle::transformProperties(T);
    U_ = transform(T, U_);
}


void Foam::bubble::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}


// ************************************************************************* //
