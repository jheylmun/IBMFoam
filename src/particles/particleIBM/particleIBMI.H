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


// * * * * * * * * * * * * * * * * Access funtions  * * * * * * * * * * * * * * //

Foam::label Foam::particleIBM::index() const
{
    return index_;
}


Foam::particleShape::locationType Foam::particleIBM::onMesh() const
{
    return shape_->onMesh();
}


bool Foam::particleIBM::active() const
{
    return active_;
}


const Foam::vector& Foam::particleIBM::center() const
{
    return shape_->center_;
}


Foam::vector& Foam::particleIBM::center()
{
    return shape_->center_;
}


Foam::scalar Foam::particleIBM::r(const vector& x2) const
{
    return shape_->r(x2);
}


Foam::scalar Foam::particleIBM::d() const
{
    return shape_->d();
}


Foam::scalar Foam::particleIBM::mass() const
{
    return shape_->V()*rho_;
}


const Foam::vector& Foam::particleIBM::v() const
{
    return v_;
}


Foam::vector& Foam::particleIBM::v()
{
    return v_;
}


const Foam::vector& Foam::particleIBM::omega() const
{
    return omega_;
}


Foam::vector& Foam::particleIBM::omega()
{
    return omega_;
}


const Foam::scalar& Foam::particleIBM::rho() const
{
    return rho_;
}


Foam::scalar& Foam::particleIBM::rho()
{
    return rho_;
}


const Foam::vector& Foam::particleIBM::F() const
{
    return integratedForce_;
}


Foam::vector& Foam::particleIBM::F()
{
    return integratedForce_;
}


const Foam::scalar& Foam::particleIBM::age() const
{
    return age_;
}

Foam::scalar& Foam::particleIBM::age()
{
    return age_;
}

Foam::particleIBM* Foam::particleIBM::particlePtr()
{
    return particlePtr_;
}

Foam::label Foam::particleIBM::centerProc() const
{
    return shape_->centerProc_;
}

const Foam::boolList& Foam::particleIBM::neiProcs() const
{
    return shape_->neiProcs_;
}

// ************************************************************************* //