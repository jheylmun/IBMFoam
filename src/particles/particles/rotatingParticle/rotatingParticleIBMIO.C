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

#include "rotatingParticleIBM.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class pType>
const std::size_t Foam::rotatingParticleIBM<pType>::sizeofFields_
(
    sizeof(rotatingParticleIBM<pType>)
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class pType>
Foam::rotatingParticleIBM<pType>::rotatingParticleIBM
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    pType(mesh, is, readFields),
    omega_(Zero)
{}

// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

template<class pType>
void Foam::rotatingParticleIBM<pType>::writeFields
(
    const typename Cloud<rotatingParticleIBM<pType>>& c
)
{
    pType::writeFields(c);
    label np = c.size();

    IOField<vector> omega(c.fieldIOobject("omega", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename Cloud<rotatingParticleIBM<pType>>, c, iter)
    {
        omega[i] = pIter().omega();
        i++;
    }

    const bool valid = np > 0;
    omega.write(valid);
}


template<class pType>
void Foam::rotatingParticleIBM<pType>::readFields
(
    typename Cloud<rotatingParticleIBM<pType>>& c
)
{
    bool valid = c.size();

    pType::readFields(c);

    IOField<vector> omega
    (
        c.fieldIOobject("omega", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, omega);

    label i = 0;
    forAllIter(typename Cloud<rotatingParticleIBM<pType>>, c, iter)
    {
        pIter().omega_ = omega[i];
        i++;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class pType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const rotatingParticleIBM<pType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const pType&>(p)
            << token::SPACE << p.omega();
    }
    else
    {
        os  << static_cast<const pType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            rotatingParticleIBM<pType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const rotatingParticleIBM<pType>&)"
    );

    return os;
}


// ************************************************************************* //