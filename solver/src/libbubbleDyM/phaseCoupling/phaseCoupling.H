/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::phaseCoupling

Description
    Class for the coupling between VoF elements and lagrangian bubbleClouds
    in the corresponding bubbleCloud.


\*---------------------------------------------------------------------------*/

#ifndef phaseCoupling_H
#define phaseCoupling_H

#include "fvCFD.H"
#include "vectorList.H"
#include "bubbleCloud.H"
#include "IOdictionary.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseCoupling Declaration
\*---------------------------------------------------------------------------*/

class phaseCoupling
{
    // Private data

        //- Reference to mesh and dictionary
        const fvMesh& mesh_;
        const IOdictionary dict_;
        
        // Reference to phase fraction field
        volScalarField& beta_;
        // Reference to velocity field
        const volVectorField& U_;

        //  Reference to phase fraction field
        volScalarField& alpha_;
        
        // Reference to phase pressure
        volScalarField& p_;
        // Reference to lagrangian bubbleCloud cloud from solver
        bubbleCloud& cloud_;

	// On/off switch
        bool VoFIDactive_;
        bool alphaEqnCouple_;
    
        bool momentumCouple_;

        // Limit for alpha which counts towards bubbleClouds
        scalar alphaLimit_;

        // Start time
        scalar startTime_;

        // Interval for phase coupling
        label nInterval_;

        scalar vlim_;

        // Marker field to indicate bubbleCloud volume ID
        volScalarField vofID_;

        scalar Rthr_;
        bool   Rcriterion_;

        label  Nle_;
        bool   NleCriterion_;      
        label  Nel_;
        bool   NelCriterion_;      

        label tNumber_;

        scalarList  pRadius_;
        vectorList  pCenter_;
        vectorList  pU_;
        labelList   pGrids_;
        labelList   closeVoF_;
        vectorList  oldMom_;
        vectorList  newMom_;
        labelList pTracks_;
        scalar cVolume_; 
        scalar cddRt_; 
        labelListList  overlap_;
        scalarList  innerVolume_;
        scalarList  interfaceVolume_;

        volScalarField effectZone_;

        volVectorField momentumSource_;

        //write lagrangian bubble radius, velocity and position to file
        autoPtr<OFstream> outputPtrR_;
        autoPtr<OFstream> outputPtrV_;
        autoPtr<OFstream> outputPtrX_;
        autoPtr<OFstream> outputPtrddRt_;


        //- Check, if neighbours of cellI belong to connected volI
        void checkvolumeConnection
        (
            const label cellI,
            const label iter
        );


public:

    // Constructors

        //- Construct from components
        phaseCoupling
        (
            const volVectorField& U,
            volScalarField& alpha,
            volScalarField& p,
            volScalarField& beta,
            bubbleCloud& cloud_
        );


    // Destructor - default

    // Member Functions

        inline bool alphaEqnCouple() const;
        inline bool momentumCouple() const;

        void initial();
        Foam::word relativePos(const vector& center, const scalar& radius, const label& celli);

        Foam::scalar gaussiankernelfunction(scalar x, scalar R);

        Foam::scalar voidfraction(scalar x, scalar R, scalar vcell, scalar vb);

        void updateBubble(bool remove = true);
        void updateBeta(volScalarField& beta);
        void updateAlphaSource(volScalarField& Sa);
        const volVectorField& momentumSource();

        void writeLagrangianBubble();

        //- Detect bubbleClouds and convert to bubbleClouds
        void update();
        
        //- transition with volume higher than limit pressure
        void updatevolume();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseCouplingI.H"


#endif

// ************************************************************************* //
