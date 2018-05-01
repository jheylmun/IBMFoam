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
#include "particleShape.H"
#include "wallFvPatch.H"
#include "randomDistribution.H"

// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * * //

template<class pType>
Foam::vector Foam::movingParticleIBM<pType>::v(const vector& pt) const
{
    return v_;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class pType>
Foam::movingParticleIBM<pType>::movingParticleIBM
(
    const polyMesh& mesh,
    const dictionary& dict,
    const label index
)
:
    pType(mesh, dict, index),
    v_(dict.template lookupOrDefault<vector>("v", Zero))
{
    if (mesh_.findCell(shape_->center_) != (-1))
    {
        active_ = 1;
    }
    else
    {
        active_ = 0;
    }

    if (active_)
    {
        track(shape_->center() - position(), 0);
    }
}


template<class pType>
Foam::movingParticleIBM<pType>::movingParticleIBM
(
    const movingParticleIBM<pType>& p,
    const label copy,
    const vector& center,
    const vector& theta,
    const vector& v,
    const vector& omega,
    const scalar& age
)
:
    pType(p, copy, center, theta, v, omega, age)
    v_(v)
{}


template<class pType>
Foam::movingParticleIBM<pType>::movingParticleIBM
(
    const movingParticleIBM<pType>& p
)
:
    pType(p)
    v_(p.v_)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * *//

template<class pType>
Foam::movingParticleIBM<pType>::~movingParticleIBM()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

template<class pType>
void Foam::movingParticleIBM<pType>::solve
(
    const scalar& dt
)
{
    if (!onMesh())
    {
        return;
    }

    vector F =
        integratedForce_
      + sum(collisionForces_)
      + sum(wallForces_);

    v_ += dt*(F/mass() + g) ;
    vector dx = dt*v_;
    shape_->center_ += dx;

    if (mesh_.findCell(shape_->center_) != (-1))
    {
        active_ = 1;
    }
    else
    {
        active_ = 0;
    }

    if (active_)
    {
        track(shape_->center() - position(), 0);
    }
    }

    collisionForces_ = Zero;
    wallForces_ = Zero;
}


template<class pType>
void Foam::movingParticleIBM<pType>::update()
{
    shape_->moveMesh(shape_->center_);
    shape_->updateCellLists();

    if (mesh_.findCell(shape_->center_) != (-1))
    {
        active_ = 1;
    }
    else
    {
        active_ = 0;
    }

    if (active_)
    {
        track(shape_->center() - position(), 0);
    }
}