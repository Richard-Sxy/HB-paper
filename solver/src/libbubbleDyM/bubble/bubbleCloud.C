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
#include "fvMesh.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "fvcCurl.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubbleCloud::bubbleCloud
(
    const fvMesh& mesh,
    volScalarField& alpha,
    volScalarField& p,
    const word& cloudName,
    bool readFields
)
:
    Cloud<bubble>(mesh, cloudName, false),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            "particleProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    rhop_(dimensionedScalar("rhop", particleProperties_).value()),
    e_(dimensionedScalar("e", particleProperties_).value()),
    mu_(dimensionedScalar("mu", particleProperties_).value()),
    rndGen_(Pstream::myProcNo()),
    alpha_(alpha),
    p_(p),
    positionsFile_(particleProperties_.subDict("distribution").lookup("positionsFile")),
    positions_
    (
        IOobject
	    (
           positionsFile_,
	       this->db().time().constant(),
	       mesh_,
	       IOobject::MUST_READ_IF_MODIFIED,
	       IOobject::NO_WRITE
	    )
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
           particleProperties_.subDict("distribution").subDict("sizeDistribution"), rndGen_
        )
    ),
    DO_(0.0), 
    RO_(0.0),
    p0_(dimensionedScalar("p0", particleProperties_).value()),
    tInjStart_(dimensionedScalar("tInjStart", particleProperties_).value()),
    tInjStartWall_(dimensionedScalar("tInjStartWall", particleProperties_).value()),
    tInjEnd_(dimensionedScalar("tInjEnd", particleProperties_).value()),
    tInjEndWall_(dimensionedScalar("tInjEndWall", particleProperties_).value()),
    injectInterval_(dimensionedScalar("injectInterval", particleProperties_).value()),
    pSat_(dimensionedScalar("pSat", particleProperties_).value()),
    ranGen_(clock::getTime() + pid()),
    deltaTL_(particleProperties_.subDict("particleProperty").lookupOrDefault("deltaTL",1e-8)),
    sigma_(dimensionedScalar("sigma", particleProperties_).value()),
    gamma_(readScalar(particleProperties_.subDict("particleProperty").lookup("gamma"))),
    solveRP_(particleProperties_.subDict("particleProperty").lookupOrDefault<Switch>("solveRP","true")),
    randomNum_(particleProperties_.subDict("particleProperty").lookupOrDefault("randomNum",0.5)),
    randomInject_(particleProperties_.subDict("particleProperty").lookupOrDefault<Switch>("randomInject","true"))
{
    if (readFields)
    {
        bubble::readFields(*this);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::bubbleCloud::move(bubbleCloud& cloud, const dimensionedVector& g)
{
    const volScalarField& rho = mesh_.lookupObject<const volScalarField>("rho");
    const volVectorField& U = mesh_.lookupObject<const volVectorField>("U");
    const volScalarField& nu = mesh_.lookupObject<const volScalarField>("nu");
    const volScalarField& p = mesh_.lookupObject<const volScalarField>("p");
    const volScalarField alphaC = mesh_.lookupObject<const volScalarField>("alpha.water");
    const volVectorField& curlUc = mesh_.lookupObject<const volVectorField>("curlUc");
    const volVectorField& DUcDt = mesh_.lookupObject<const volVectorField>("DUcDt"); 
    const volVectorField& GradPc = mesh_.lookupObject<const volVectorField>("GradPc");   

    interpolationCellPoint<scalar> rhoInterp(rho);
    interpolationCellPoint<vector> UInterp(U);
    interpolationCellPoint<scalar> nuInterp(nu);
    interpolationCellPoint<scalar> pInterp(p);
    interpolationCellPoint<scalar> alphaCInterp(alphaC);
    interpolationCellPoint<vector> curlUcInterp(curlUc);
    interpolationCellPoint<vector> DUcDtInterp(DUcDt);
    interpolationCellPoint<vector> GradPcInterp(GradPc);

    // inject in area that the pressure is smaller than pSat
    bubble::trackingData
        td(*this, rhoInterp, UInterp, nuInterp, pInterp,alphaCInterp,curlUcInterp,DUcDtInterp, GradPcInterp,g.value());

    dimensionedScalar deltaT = mesh_.time().deltaT();
    scalar TL = 0.0;
    scalar deltaTL = deltaTL_;
    if(deltaT.value() < deltaTL_)
    {
        Cloud<bubble>::move(*this, td, mesh_.time().deltaTValue());
    if(Pstream::master())
    {
        Info<<"bubble deltaT  : "<<mesh_.time().deltaTValue()<<endl;
    }
    }
    else 
    {
        while((deltaT.value() - TL) > 1e-8)
        {
            if((deltaT.value()-TL)<deltaTL)
            {
                deltaTL = deltaT.value()-TL;
            }
            Cloud<bubble>::move(*this, td, deltaTL);
            TL += deltaTL;
        }
    }

    if(Pstream::master())
    {
        Info<<"Update bubble motion and Radius successful"<<endl;
    }
    // inject from mean flow 
    this->inject(td);

    // inject the bubble in the surface of solid gap
    this->injectwall(td);

    this->remove();
    label nTotParcel = returnReduce(this->size(), sumOp<label>());
    Info << "BubbleCloud size  = " << nTotParcel<<endl;
}

void Foam::bubbleCloud::inject
(
    bubble::trackingData &td
)
{
    if(mesh_.time().value()> tInjStart_ && mesh_.time().value()< tInjEnd_)
    {
        forAll(positions_,pI)
        {
            scalar x = 0.0;
            if(randomInject_)
            {
               x = ranGen_.sample01<scalar>();
            }
            if (x < randomNum_)
            {
                   label cellI = -1;
                  label tetFaceI = -1;
                    label tetPtI = -1;
	          vector posB = positions_[pI];
                 vector velocity = vector::zero; 
	          mesh_.findCellFacePt(posB, cellI, tetFaceI, tetPtI);
                 label proci = -1.0;
/*            
                  if (cellI >= 0)
                 { 
       		         proci = Pstream::myProcNo();
    	         }
            
    	           reduce(proci, maxOp<label>());

            // Ensure that only one processor attempts to insert this Parcel
            if (proci != Pstream::myProcNo())
   	    {
       		 cellI = -1;
       		 tetFaceI = -1;
      		  tetPtI = -1;
    	    }
           // Last chance - find nearest cell and try that one - the point is
           //      probably on an edge
                if (proci == -1)
         {
            cellI = mesh_.findNearestCell(posB);

           if (cellI >= 0)
            {
                posB += 1e-3*(mesh_.C()[cellI] - posB);

              mesh_.findCellFacePt
            (
                posB,
                cellI,
                tetFaceI,
                tetPtI
            );

            if (cellI > 0)
            {
                proci = Pstream::myProcNo();
            }
        }

        reduce(proci, maxOp<label>());

        if (proci != Pstream::myProcNo())
        {
            cellI = -1;
            tetFaceI = -1;
            tetPtI = -1;
        }
    }
 */       //-----    出现重复的区域可能位于Proc边界 因此可以不考虑目前
	        if(cellI > 0)
            {    
                 scalar diameter = sizeDistribution_().sample();
                 bubble* ptr1 = new bubble(mesh_, posB,cellI, tetFaceI, tetPtI, diameter, velocity);
                 Cloud<bubble>::addParticle(ptr1);
            }
	}
       } 
    tInjStart_ += injectInterval_;
    Info<<"Inject bubble from flow successful "<<endl;
  } 
}


void Foam::bubbleCloud::injectwall
(
    bubble::trackingData &td
)
{
     // 更合适的方式
    vector wallPos1_(particleProperties_.subDict("particleProperty").lookup("wallPosition1"));
    vector wallPos2_(particleProperties_.subDict("particleProperty").lookup("wallPosition2"));
    if(mesh_.time().value() > tInjStartWall_  && mesh_.time().value() < tInjEndWall_)
    {
          forAll(mesh_.cells(),cellI)
          {
               if (alpha_[cellI] > 0.6 && p_[cellI]  < 2000.0)
               {
                  vector position = mesh_.cellCentres()[cellI];
                  scalar x = ranGen_.sample01<scalar>();
                  if((wallPos1_.x() <position[0] && position[0]< wallPos2_.x()) && (wallPos1_.y()<position[1] && position[1]<wallPos2_.y()))
	              { 
	                if( x <= 0.5)
	                {
	    	          scalar diameter = sizeDistribution_().sample();
	   	              vector velocity = vector::zero;
		              label celli;
                      label tetFacei;
                      label tetPti;

                      mesh_.findCellFacePt(position,celli,tetFacei, tetPti);

	                  if(celli > 0)
                      {
                           bubble* ptr = new bubble(mesh_, position,celli, tetFacei, tetPti, diameter, velocity);
	                       Cloud<bubble>::addParticle(ptr);
                      }
                    }
                  }     
              }
          }
          tInjStartWall_ += injectInterval_;
    Info<<"Inject bubble from Wall successful "<<endl;
    }
}
void Foam::bubbleCloud::injectbubble
(
     vector position,
     const scalar diameter,
     const vector velocity
)
{
     label cellI;
     label tetFaceI;
     label tetPtI;

     mesh_.findCellFacePt(position,cellI,tetFaceI, tetPtI);
     label proci = -1;

            if (cellI >= 0)
            {
                 proci = Pstream::myProcNo();
            }

     reduce(proci, maxOp<label>());
     
             if (proci != Pstream::myProcNo())
            {
                 cellI = -1;
                 tetFaceI = -1;
                  tetPtI = -1;
            }  
            if (proci == -1)
         {
            cellI = mesh_.findNearestCell(position);

           if (cellI >= 0)
            {
                position += 1e-3*(mesh_.C()[cellI] - position);

              mesh_.findCellFacePt
            (
                position,
                cellI,
                tetFaceI,
                tetPtI
            );

            if (cellI > 0)
            {
                proci = Pstream::myProcNo();
            }
        } 
} 
     scalar R0 = 1e-6;
	 if(cellI > 0)
     {
         bubble* ptr = new bubble(mesh_, position,cellI, diameter, R0, velocity);
	     Cloud<bubble>::addParticle(ptr);
     }
}
void Foam::bubbleCloud::remove()
{
    vector failureArea_(particleProperties_.subDict("particleProperty").lookup("failureArea"));
    scalar injectAreaVolume_(readScalar(particleProperties_.subDict("particleProperty").lookup("injectAreaVolume")));   
    DO_ = 0.0;
    RO_ = 0.0;
    for(bubble& b : *this)
    {
        DO_ += 4.0/3.0*constant::mathematical::pi*pow(b.R0(),3.0)/(injectAreaVolume_+VSMALL);
        RO_ += 4.0/3.0*constant::mathematical::pi*pow(b.R(),3.0)/(injectAreaVolume_+VSMALL);
        if(mag(b.position().x()) > mag(failureArea_.x()))
        {
            DO_ -= 4.0/3.0*constant::mathematical::pi*pow(b.R0(),3.0)/(injectAreaVolume_+VSMALL);
            RO_ -= 4.0/3.0*constant::mathematical::pi*pow(b.R(),3.0)/(injectAreaVolume_+VSMALL);
            this->deleteParticle(b);
        }
    }
    reduce(DO_,sumOp<scalar>());
    reduce(RO_,sumOp<scalar>());
    if(Pstream::master())
    {
        Info<<"RO: "<<RO_<<" DO: "<<DO_<<endl;
    }
}
// ************************************************************************* //
// ************************************************************************* //
