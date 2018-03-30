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

#include "particleCloud.H"
#include "particleIBM.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::particleIBM::sizeofFields_
(
    sizeof(particleIBM)
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleIBM::particleIBM
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    mesh_(mesh),
    dict_(mesh.lookupObject<IOdictionary>("initialStates").subDict("default")),
    index_(-1),
    active_(0),
    shape_(particleShape::New(mesh_, dict_, position())),
    v_(Zero),
    omega_(Zero),
    rho_(0.0),
    age_(0.0),
    integratedForce_(Zero)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            active_ = readBool(is);
            is >> v_;
            is >> omega_;
            rho_ = readScalar(is);
            age_ = readScalar(is);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&active_), sizeofFields_);
        }
    }

    // Check state of Istream
    is.check
    (
        "particleIBM::patricleIBM"
        "(const polyMesh&, const dictinoary&, Istream&, bool)"
    );
}

// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::particleIBM::writeFields(const Cloud<particleIBM>& c)
{
    particle::writeFields(c);
    label np = c.size();

    IOField<label> index(c.fieldIOobject("index", IOobject::NO_READ), np);
    IOField<label> active(c.fieldIOobject("active", IOobject::NO_READ), np);
    IOField<vector> center(c.fieldIOobject("center", IOobject::NO_READ), np);
    IOField<label> centerI(c.fieldIOobject("centerI", IOobject::NO_READ), np);
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<vector> v(c.fieldIOobject("v", IOobject::NO_READ), np);
    IOField<vector> omega(c.fieldIOobject("omega", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);
    IOField<scalar> age(c.fieldIOobject("age", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename Cloud<particleIBM>, c, iter)
    {
        const particleIBM& p = iter();

        index[i] = p.index();
        active[i] = p.active();
        center[i] = p.center();
        centerI[i] = p.centerIndex();
        d[i] = p.d();
        v[i] = p.v();
        omega[i] = p.omega();
        rho[i] = p.rho();
        age[i] = p.age();
        i++;
    }

    const bool valid = np > 0;

    index.write(valid);
    active.write(valid);
    center.write(valid);
    centerI.write(valid);
    d.write(valid);
    v.write(valid);
    omega.write(valid);
    rho.write(valid);
    age.write(valid);
}


void Foam::particleIBM::readFields(Cloud<particleIBM>& c)
{
    bool valid = c.size();

    particle::readFields(c);

    IOField<label> index
    (
        c.fieldIOobject("index", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, index);
    IOField<label> active
    (
        c.fieldIOobject("active", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, active);

    IOField<scalar> d
    (
        c.fieldIOobject("d", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, d);

    IOField<scalar> theta
    (
        c.fieldIOobject("theta", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, theta);

    IOField<vector> v
    (
        c.fieldIOobject("v", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, v);

    IOField<vector> omega
    (
        c.fieldIOobject("omega", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, omega);

    IOField<scalar> rho
    (
        c.fieldIOobject("rho", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, rho);

    IOField<scalar> age
    (
        c.fieldIOobject("age", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, age);


    label i = 0;

    forAllIter(typename Cloud<particleIBM>, c, iter)
    {
        particleIBM& p = iter();

        p.index_ = index[i];
        p.active_ = active[i];
        p.v_ = v[i];
        p.omega_ = omega[i];
        p.rho_ = rho[i];
        p.age_ = age[i];

        i++;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const particleIBM& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.index()
            << token::SPACE << p.active()
            << token::SPACE << p.position()
            << token::SPACE << p.v()
            << token::SPACE << p.omega()
            << token::SPACE << p.rho()
            << token::SPACE << p.writeShape(os);
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            particleIBM::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const particleIBM&)"
    );

    return os;
}


Foam::Ostream& Foam::particleIBM::writeShape(Ostream& os) const
{
    os << shape_();
    return os;
}

// ************************************************************************* //