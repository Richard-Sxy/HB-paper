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

#include "bubble.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::bubble::sizeofFields
(
    sizeof(bubble) - sizeof(particle)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubble::bubble
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is  >> d_ >> U_ >> R_>>dRt_>>ddRt_>>R0_>>U0_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, &d_);
            readRawScalar(is, U_.data(), vector::nComponents);
            readRawScalar(is, U0_.data(), vector::nComponents);
            readRawScalar(is, &R_);
            readRawScalar(is, &R0_);
            readRawScalar(is, &dRt_);
            readRawScalar(is, &ddRt_);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&d_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


void Foam::bubble::readFields(Cloud<bubble>& c)
{
    bool valid = c.size();

    particle::readFields(c);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, d);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ), valid);
    IOField<vector> U0(c.fieldIOobject("U0", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, U);
    c.checkFieldIOobject(c, U0);

    IOField<scalar> R(c.fieldIOobject("R", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, R);
    IOField<scalar> R0(c.fieldIOobject("R0", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, R0);
    IOField<scalar> dRt(c.fieldIOobject("dRt", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, dRt);
    IOField<scalar> ddRt(c.fieldIOobject("ddRt", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, ddRt);
    
    label i = 0;
    for (bubble& p : c)
    {
        p.d_ = d[i];
        p.U_ = U[i];
        p.R_ = R[i];
        p.U0_ = U0[i];
        p.R0_ = R0[i];
        p.dRt_ = dRt[i];
        p.ddRt_ = ddRt[i];
        ++i;
    }
}


void Foam::bubble::writeFields(const Cloud<bubble>& c)
{
    particle::writeFields(c);

    const label np = c.size();

    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<vector> U0(c.fieldIOobject("U0", IOobject::NO_READ), np);
    IOField<scalar> R(c.fieldIOobject("R", IOobject::NO_READ), np);
    IOField<scalar> R0(c.fieldIOobject("R0", IOobject::NO_READ), np);
    IOField<scalar> dRt(c.fieldIOobject("dRt", IOobject::NO_READ), np);
    IOField<scalar> ddRt(c.fieldIOobject("ddRt", IOobject::NO_READ), np);
    label i = 0;
    for (const bubble& p : c)
    {
        d[i] = p.d_;
        U[i] = p.U_;
        R[i] = p.R_;
        U0[i] = p.U0_;
        R0[i] = p.R0_;
        dRt[i] = p.dRt_;
        ddRt[i] = p.ddRt_;
        ++i;
    }

    d.write(np > 0);
    U.write(np > 0);
    R.write(np > 0);
    R0.write(np > 0);
    U0.write(np > 0);
    dRt.write(np > 0);
    ddRt.write(np > 0);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const bubble& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.d_
            << token::SPACE << p.U_
            << token::SPACE << p.U0_
            << token::SPACE << p.R0_
            << token::SPACE << p.dRt_
            << token::SPACE << p.ddRt_
            << token::SPACE << p.R_;

    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.d_),
            bubble::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
