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

#include "nonRotatingParticleIBM.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class pType>
const std::size_t Foam::nonRotatingParticleIBM<pType>::sizeofFields_
(
    sizeof(nonRotatingParticleIBM<pType>)
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class pType>
Foam::nonRotatingParticleIBM<pType>::nonRotatingParticleIBM
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    pType(mesh, is, readFields)
{}

// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

template<class pType>
void Foam::nonRotatingParticleIBM<pType>::writeFields
(
    const typename Cloud<nonRotatingParticleIBM<pType>>& c
)
{
    pType::writeFields(c);
}


template<class pType>
void Foam::nonRotatingParticleIBM<pType>::readFields
(
    typename Cloud<nonRotatingParticleIBM<pType>>& c
)
{
    pType::readFields(c);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class pType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const nonRotatingParticleIBM<pType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const pType&>(p);
    }
    else
    {
        os  << static_cast<const pType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            nonRotatingParticleIBM<pType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const nonRotatingParticleIBM&)"
    );

    return os;
}

// ************************************************************************* //