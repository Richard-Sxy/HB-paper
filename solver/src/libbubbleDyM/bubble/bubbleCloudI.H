/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::bubbleCloud::mesh() const
{
    return mesh_;
}


inline Foam::scalar Foam::bubbleCloud::rhop() const
{
    return rhop_;
}


inline Foam::scalar Foam::bubbleCloud::e() const
{
    return e_;
}
inline Foam::scalar Foam::bubbleCloud::mu() const
{
    return mu_;
}


inline Foam::scalar Foam::bubbleCloud::pSat() const
{
return pSat_;
}

inline Foam::scalar Foam::bubbleCloud::sigma() const
{
return sigma_;
}
inline Foam::Random& Foam::bubbleCloud::rndGen() const
{
return rndGen_;
}
inline Foam::scalar Foam::bubbleCloud::p0() const
{
return p0_;
}
inline bool  Foam::bubbleCloud::solveRP() const
{
    return solveRP_;
}
inline Foam::scalar Foam::bubbleCloud::gamma() const
{
    return gamma_;
}

inline Foam::scalar Foam::bubbleCloud::deltaTL() const
{
    return deltaTL_;
}
// ************************************************************************* //
