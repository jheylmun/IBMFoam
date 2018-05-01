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

#include "stationaryParticleIBM.H"
#include "particleShape.H"
#include "wallFvPatch.H"
#include "randomDistribution.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class pType>
Foam::vector Foam::stationaryParticleIBM<pType>::v(const vector& pt) const
{
    return Zero;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class pType>
Foam::stationaryParticleIBM<pType>::stationaryParticleIBM
(
    const polyMesh& mesh,
    const dictionary& dict,
    const label index
)
:
    pType(mesh, dict, index)
{}


template<class pType>
Foam::stationaryParticleIBM<pType>::stationaryParticleIBM
(
    const stationaryParticleIBM<pType>& p,
    const label copy,
    const vector& center,
    const vector& theta,
    const vector& v,
    const vector& omega,
    const scalar& age
)
:
    pType(p, copy, center, theta, v, omega, age)
{}


template<class pType>
Foam::stationaryParticleIBM<pType>::stationaryParticleIBM
(
    const stationaryParticleIBM<pType>& p
)
:
    pType(p)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * *//

template<class pType>
Foam::stationaryParticleIBM::~stationaryParticleIBM()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

template<class pType>
void Foam::stationaryParticleIBM<pType>::solve
(
    const scalar& dt,
    const vector& g,
    const bool moving,
    const bool rotation
)
{
    return;
}


void Foam::stationaryParticleIBM<pType>::update()
{
    return;
}