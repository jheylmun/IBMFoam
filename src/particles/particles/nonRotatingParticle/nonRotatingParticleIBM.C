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

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class pType>
Foam::vector Foam::nonRotatingParticleIBM<pType>::v(const vector& pt) const
{
    return (-(pt - shape_->center_)^omega_ ) + v();
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class pType>
Foam::nonRotatingParticleIBM<pType>::nonRotatingParticleIBM
(
    const polyMesh& mesh,
    const dictionary& dict,
    const label index
)
:
    pType(mesh, dict, index)
{}


template<class pType>
Foam::nonRotatingParticleIBM<pType>::nonRotatingParticleIBM
(
    const nonRotatingParticleIBM<pType>& p,
    const label copy,
    const vector& center,
    const vector& theta,
    const vector& v,
    const vector& omega,
    const scalar& age
)
:
    pType(p,copy, center, theta, v, omega, age)
{}


template<class pType>
Foam::nonRotatingParticleIBM<pType>::nonRotatingParticleIBM
(
    const nonRotatingParticleIBM<pType>& p
)
:
    pType(p)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * *//

template<class pType>
Foam::nonRotatingParticleIBM::~nonRotatingParticleIBM()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

template<class pType>
void Foam::nonRotatingParticleIBM<pType>::solve
(
    const scalar& dt
)
{
    return;
}


template<class pType>
void Foam::nonRotatingParticleIBM::update()
{}