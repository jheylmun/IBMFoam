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

#include "isothermalParticleIBM.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::isothermalParticleIBM<pType>::sizeofFields_
(
    sizeof(isothermalParticleIBM)
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isothermalParticleIBM<pType>::isothermalParticleIBM
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    pType(mesh, is, readFields),
    mesh_(mesh),
    index_(-1),
    active_(0),
    shape_(),
    v_(Zero),
    omega_(Zero),
    rho_(0.0),
    age_(0.0),
    integratedForce_(Zero),
    integratedTorque_(Zero)
{}

// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::isothermalParticleIBM<pType>::writeFields
(
    const typename Cloud<isothermalParticleIBM<pType>>& c
)
{
    pType::writeFields(c);
    label np = c.size();

    IOField<scalar> T(c.fieldIOobject("T", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename Cloud<isothermalParticleIBM<pType>>, c, iter)
    {
        T[i] = pIter().T();
        i++;
    }

    const bool valid = np > 0;
    T.write(valid);
}


void Foam::isothermalParticleIBM<pType>::readFields
(
    typename Cloud<isothermalParticleIBM<pType>>& c
)
{
    bool valid = c.size();

    pType::readFields(c);

    IOField<scalar> T
    (
        c.fieldIOobject("T", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, T);

    label i = 0;
    forAllIter(typename Cloud<isothermalParticleIBM<pType>>, c, iter)
    {
        pIter().T_ = T[i];
        i++;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const isothermalParticleIBM<pType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const pType&>(p)
            << token::SPACE << p.T();
    }
    else
    {
        os  << static_cast<const pType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            isothermalParticleIBM<pType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const isothermalParticleIBM<pType>&)"
    );

    return os;
}

// ************************************************************************* //