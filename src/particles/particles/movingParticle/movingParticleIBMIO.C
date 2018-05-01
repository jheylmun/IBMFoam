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

#include "movingParticleIBM.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class pType>
const std::size_t Foam::movingParticleIBM<pType>::sizeofFields_
(
    sizeof(movingParticleIBM<pType>)
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class pType>
Foam::movingParticleIBM<pType>::movingParticleIBM
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    pType(mesh, is, readFields),
    v_(Zero)
{}

// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

template<class pType>
void Foam::movingParticleIBM<pType>::writeFields
(
    const typename Cloud<movingParticleIBM<pType>>& c
)
{
    pType::writeFields(c);
    label np = c.size();

    IOField<vector> v(c.fieldIOobject("v", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename Cloud<movingParticleIBM>, c, iter)
    {
        v[i] = pIter.v();
        i++;
    }

    const bool valid = np > 0;

    v.write(valid);
}


template<class pType>
void Foam::movingParticleIBM<pType>::readFields
(
    typename Cloud<movingParticleIBM<pType>>& c
)
{
    bool valid = c.size();

    pType::readFields(c);

    IOField<vector> v
    (
        c.fieldIOobject("v", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, v);

    const polyMesh& mesh = c.pMesh();

    label i = 0;
    forAllIter(typename Cloud<movingParticleIBM<pType>>, c, iter)
    {
        pIter().v_ = v[i];
        i++;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class pType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const movingParticleIBM<pType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const typename pType&>(p)
            << p.v();
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            movingParticleIBM<pType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const movingParticleIBM<pType>&)"
    );

    return os;
}


// ************************************************************************* //