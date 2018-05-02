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

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class pType>
Foam::vector Foam::rotatingParticleIBM<pType>::v(const vector& pt) const
{
    return ((shape_->center_ - pt)^omega_ ) + pType::v(pt);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class pType>
Foam::rotatingParticleIBM<pType>::rotatingParticleIBM
(
    const polyMesh& mesh,
    const dictionary& dict,
    const label index
)
:
    pType(mesh, dict, index),
    omega_(dict.template lookupOrDefault<vector>("omega", Zero))
{}


template<class pType>
Foam::rotatingParticleIBM<pType>::rotatingParticleIBM
(
    const rotatingParticleIBM<pType>& p,
    const label copy,
    const vector& center,
    const vector& theta,
    const vector& v,
    const vector& omega,
    const scalar& age
)
:
    pType(p,copy, center, theta, v, omega, age),
    omega_(omega)
{}


template<class pType>
Foam::rotatingParticleIBM<pType>::rotatingParticleIBM
(
    const rotatingParticleIBM<pType>& p
)
:
    pType(p),
    omega_(p.omega_)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * *//

template<class pType>
Foam::rotatingParticleIBM::~rotatingParticleIBM()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

template<class pType>
void Foam::rotatingParticleIBM<pType>::solve
(
    const scalar& dt
)
{
    pType::solve(dt);

    vector T =
        integratedTorque_
      + sum(collisionTorques_)
      + sum(wallTorques_);

    omega_ += dt*T/(mass()*shape_->I());
    shape_->theta() += omega_*dt;

    collisionTorques_ = Zero;
    wallTorques_ = Zero;
}


template<class pType>
void Foam::rotatingParticleIBM::update()
{
    pType::update();
}