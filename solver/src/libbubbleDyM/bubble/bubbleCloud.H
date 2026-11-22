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

Class
    Foam::bubbleCloud

Description
    A Cloud of bubble

SourceFiles
    bubbleCloudI.H
    bubbleCloud.C
    bubbleCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef bubbleCloud_H
#define bubbleCloud_H

#include "Cloud.H"
#include "bubble.H"
#include "IOdictionary.H"
#include "Random.H"
#include "globalIOFields.H"
#include "distributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class bubbleCloud Declaration
\*---------------------------------------------------------------------------*/

class bubbleCloud
:
    public Cloud<bubble>
{
    // Private data

        const fvMesh& mesh_;

        IOdictionary particleProperties_;

        scalar rhop_;
        scalar e_;
	    scalar mu_;
        
        Random ranGen_;
        //- Parcel size distribution model
        const autoPtr<distributionModel> sizeDistribution_;
 
        //- Random number generator - used by injection
        mutable Random rndGen_;

        volScalarField& p_;       
        volScalarField& alpha_;       
        scalar p0_;
        scalar DO_;
        scalar RO_;
        scalar tInjStart_;
        scalar tInjStartWall_;
        scalar tInjEnd_;
        scalar tInjEndWall_;
        const scalar injectInterval_;
        scalar deltaTL_;	
        scalar pSat_; // saturation pressure
        scalar sigma_; // surface tenstion
        bool   solveRP_;
        bool   randomInject_;
        scalar   randomNum_;
        scalar gamma_;
	    //- Name of file containing positions data
	    const word positionsFile_;
        
	    //- Filed of parcel positions
	    vectorGlobalIOField positions_;


    // Private Member Functions

        //- No copy construct
        bubbleCloud(const bubbleCloud&) = delete;

        //- No copy assignment
        void operator=(const bubbleCloud&) = delete;



public:

    // Constructors

        //- Construct given mesh
        explicit bubbleCloud
        (
            const fvMesh& mesh,
            volScalarField& alpha,
            volScalarField& p,
            const word& cloudName = cloud::defaultName,
            bool readFields = true
        );


    // Member Functions

        // Access

            inline const fvMesh& mesh() const;

            inline scalar rhop() const;
            inline scalar e() const;
            inline scalar mu() const;

	        inline scalar p0() const;
	        inline scalar tInjStart() const;
            inline scalar tInjEnd() const;
            inline scalar tInjEndWall() const;
            inline scalar tInjStartWall() const;
	        inline scalar pSat() const;
            inline scalar sigma() const;
            inline bool solveRP() const;
            inline scalar gamma() const;
            inline scalar deltaTL() const;
            //- Return reference to the random object
            inline Random& rndGen() const;

    // Edit

            //- Move the particles under the influence of the given
            //  gravitational acceleration
            void move(bubbleCloud& ,const dimensionedVector& g);

	    //- Inject particles according to the dictionnary particleProperties
            void inject(bubble::trackingData &td);
            void injectwall(bubble::trackingData &td);
            void injectbubble(vector position,const scalar diameter,const vector velocity);
            void remove();

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bubbleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
