/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::bubble

Description
    Simple solid spherical particle class with one-way coupling with the
    continuous phase.

SourceFiles
    bubbleI.H
    bubble.C
    bubbleIO.C

\*---------------------------------------------------------------------------*/

#ifndef bubble_H
#define bubble_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"
#include "mathematicalConstants.H"
#include "ODESystem.H"
#include "ODESolver.H"

using namespace Foam::constant::mathematical;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class bubbleCloud;
class bubble;

Ostream& operator<<(Ostream&, const bubble&);


/*---------------------------------------------------------------------------*\
                           Class bubble Declaration
\*---------------------------------------------------------------------------*/

class bubble
:
    public particle
{
    // Private Data

        //- Diameter
        scalar d_;

        //- Velocity of parcel
        vector U_;
        vector U0_;
	    scalar R_; // -radius
        scalar dRt_; // -r_dot
        scalar R0_; // initial equilibrium radius
        scalar ddRt_;
        // solver ode 
        class RPODE  
        :
             public ODESystem
        {
                scalar p_;
                scalar p0_;
                scalar rho_;
		        scalar mu_;
		        scalar pv_;
		        scalar sigma_;
                scalar R0_;
        public:
                // Constructors
                inline RPODE
                (
                    scalar p,
                    scalar p0,
                    scalar rho,
 		            scalar mu,
		            scalar pSat,
		            scalar sigma,
                    scalar R0
                );
             
                // Member Functions
                //- Return the number of equations in the system
                inline label nEqns() const;
 
                //- Calculate the derivatives in dydx
                inline void derivatives
                (
                     const scalar x,
                     const scalarField& y,
                     scalarField& dydx
                ) const;

                //- Calculate the Jacobian of the system
                //  Need by the stiff-system solvers
                inline void jacobian
                (
                    const scalar x,
                    const scalarField& y,
                    scalarField& dfdx,                                                                                                                     
                    scalarSquareMatrix& dfdy
                ) const;
 
        };

public:

    friend class Cloud<bubble>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& rhoInterp_;
            const interpolationCellPoint<vector>& UInterp_;
            const interpolationCellPoint<scalar>& nuInterp_;
	        const interpolationCellPoint<scalar>& pInterp_;
	        const interpolationCellPoint<scalar>& alphaCInterp_;
            const interpolationCellPoint<vector>& curlUcInterp_;
            const interpolationCellPoint<vector>& DUcDtInterp_;
	        const interpolationCellPoint<vector>& GradPcInterp_;
	    //- Local gravitational or other body-force acceleration
            const vector& g_;


    public:

        // Constructors

            inline trackingData
            (
                const bubbleCloud& spc,
                const interpolationCellPoint<scalar>& rhoInterp,
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& nuInterp,
		        const interpolationCellPoint<scalar>& pInterp,
                const interpolationCellPoint<scalar>& alphaCInterp,
		        const interpolationCellPoint<vector>& curlUcInterp,
		        const interpolationCellPoint<vector>& DUcDtInterp,
		        const interpolationCellPoint<vector>& GradPcInterp_,
		        const vector& g
            );


        // Member functions

            inline const interpolationCellPoint<scalar>& rhoInterp() const;

            inline const interpolationCellPoint<vector>& UInterp() const;

            inline const interpolationCellPoint<scalar>& nuInterp() const;

	        inline const interpolationCellPoint<scalar>& pInterp() const;

	        inline const interpolationCellPoint<scalar>& alphaCInterp() const;

	        inline const interpolationCellPoint<vector>& curlUcInterp() const;

	        inline const interpolationCellPoint<vector>& DUcDtInterp() const;

	        inline const interpolationCellPoint<vector>& GradPcInterp() const;

            inline const vector& g() const;
    };


    // Static data members

        //- Size in bytes of the fields
        static const std::size_t sizeofFields;


    // Constructors

        //- Construct from a position and a cell
        //  Searches for the rest of the required topology.
        //  Other properties are zero initialised.
        inline bubble
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli = -1
        );

        //- Construct from components

	    inline bubble
	    (
            const polyMesh& mesh,
	        const vector& position,
	        const label celli,
            const label tetFacei,
            const label tetPti,
	        const scalar d,
	        const vector& U
	    );
        inline bubble
	    (
            const polyMesh& mesh,
	        const vector& position,
	        const label celli,
	        const scalar d,
            const scalar R0,
	        const vector& U
	    );
        //- Construct from Istream
        bubble
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new bubble(*this));
        }

        //- Factory class to read-construct particles (for parallel transfer)
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<bubble> operator()(Istream& is) const
            {
                return autoPtr<bubble>::New(mesh_, is, true);
            }
        };


    // Member Functions

        // Access

            //- Return diameter
            inline scalar& d();

            //- Return velocity
            inline vector& U();
            inline vector& U0();
	        inline scalar R();

	        inline scalar R0();

	        inline scalar dRt();
            inline scalar ddRt();   

            

        // Tracking


        //- Move
        bool move(bubbleCloud&, trackingData&, const scalar);


	    scalar CdRe
            (
               scalar Re
	    ); 

	    scalar Cl
	    (
               scalar Re,
               scalar rhoc,
	           scalar muc,
               vector& curlUc
	    );
      
      // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch(bubbleCloud& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch(bubbleCloud& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(bubbleCloud& cloud, trackingData& td);

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            virtual void transformProperties(const vector& separation);


    // I-O

        static void readFields(Cloud<bubble>& c);

        static void writeFields(const Cloud<bubble>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const bubble&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for bubble
template<> struct is_contiguous<bubble> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bubbleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
