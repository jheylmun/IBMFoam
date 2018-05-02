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

#include "anisothermalParticleIBM.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class pType>
Foam::anisothermalParticleIBM<pType>::anisothermalParticleIBM
(
    const polyMesh& mesh,
    const dictionary& dict,
    const label index
)
:
    pType(mesh, dict, index),
    T_(dict.lookup("T")),
    integratedQ_(Zero)
{}


template<class pType>
Foam::anisothermalParticleIBM<pType>::anisothermalParticleIBM
(
    const anisothermalParticleIBM<pType>& p,
    const label copy,
    const vector& center,
    const vector& theta,
    const vector& v,
    const vector& omega,
    const scalar& age
)
:
    pType(p, copy, center, theta, v, omega, age),
    T_(p.T_),
    integratedQ_(p.integratedQ_)
{}


template<class pType>
Foam::anisothermalParticleIBM<pType>::anisothermalParticleIBM
(
    const anisothermalParticleIBM<pType>& p
)
:
    pType(p),
    T_(p.T_),
    integratedQ_(p.integratedQ_)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * *//

template<class pType>
Foam::anisothermalParticleIBM<pType>::~anisothermalParticleIBM()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

template<class pType>
void Foam::anisothermalParticleIBM<pType>::solve(const scalar& dt)
{
    pType::solve(dt);
}


template<class pType>
void Foam::anisothermalParticleIBM<pType>::heatFlux
(
    const surfaceVectorField& Tf,
    const surfaceVectorField& Tfold
    volVectorField& Q
) const
{
// Calcuate heat flux from particle to fluid
//     label N = shape_->N();
//
//     vectorField interpolatedT(N, Zero);
//     vectorField interpolatedTold(N, Zero);
//     vectorField interpolatedT(N, Zero);
//
//     interpolateFromMesh(Tf,interpolatedT);
//     interpolateFromMesh(Tfold,interpolatedTold);
//
//     interpolateQToMesh(interpolatedT, interpolatedTold, Q);
}


// template<class pType>
// void Foam::anisothermalParticleIBM<pType>::integrateSurfaceHeatFlux
// (
//     const surfaceSymmTensorField& tauf,
//     const surfaceScalarField& pf
// )
// {
//     symmTensorField interpolatedTau(shape_->N(), Zero);
//     scalarField interpolatedP(shape_->N(), Zero);
//
//     interpolateFromMesh(tauf,interpolatedTau);
//     interpolateFromMesh(pf,interpolatedP);
//     tmp<vectorField> tmpSf(shape_->Sf());
//     const vectorField& Sf = tmpSf();
//
//     integratedForce_ = Zero;
//     integratedTorque_ = Zero;
//
//     for (label i = 0; i < shape_->nTheta(); i++)
//     {
//         for (label k = 0; k < shape_->nk(); k++)
//         {
//             vector R =
//                 shape_->baseMesh_[shape_->index2(i,k)] - shape_->center_;
//
//             const vector& Sfi = Sf[shape_->index2(i,k)];
//             scalar magSf = mag(Sfi);
//             vector norm = Sfi/magSf;
//
//             const symmTensor& tau = interpolatedTau[shape_->index(1,i,k)];
//             vector F =
//                 interpolatedP[shape_->index(1,i,k)]
//                *shape_->Sf_[shape_->index2(i,k)]
//               + (tau & Sfi);
//             vector Fn = (F & norm)*norm;
//             vector Ft = F - Fn;
//
//             integratedForce_ += F;
//             integratedTorque_ -= Ft^R;
//         }
//     }
// }


template<class pType>
void Foam::anisothermalParticleIBM<pType>::update()
{
    pType::update();
}