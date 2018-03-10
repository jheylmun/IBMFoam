/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "particleList.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::IBM::particle::sizeofFields_
(
    sizeof(particle)
);

// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::IBM::particle::writeFields(const particleList& c)
{
    const PtrList<particle>& pList = c;
    label np = pList.size();

    IOField<label> active(c.fieldIOobject("active", IOobject::NO_READ), np);
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<vector> v(c.fieldIOobject("v", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(PtrList<IBM::particle>, pList, iter)
    {
        const particle& p = iter();

        active[i] = p.active();
        d[i] = p.d();
        v[i] = p.v();
        rho[i] = p.rho();
        i++;
    }

    const bool valid = np > 0;

    active.write(valid);
    d.write(valid);
    v.write(valid);
    rho.write(valid);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::IBM::operator<<
(
    Ostream& os,
    const IBM::particle& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << p.position()
            << token::SPACE << p.active()
            << token::SPACE << p.d()
            << token::SPACE << p.v()
            << token::SPACE << p.rho();
    }
    else
    {
        os  << static_cast<const IBM::particle&>(p);
        os.write
        (
            "d",
            IBM::particle::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const particle&)"
    );

    return os;
}


// ************************************************************************* //