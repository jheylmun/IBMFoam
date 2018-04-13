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

#include "particleIBM.H"
#include "particleShape.H"
#include "wallFvPatch.H"
#include "randomDistribution.H"

namespace Foam
{
    defineTypeNameAndDebug(particleIBM, 0);
    defineTemplateTypeNameAndDebug(Cloud<particleIBM>, 0);
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::particleIBM::interpolateFToMesh
(
    const vectorList& U,
    const vectorList& Uold,
    const vectorList& S,
    volVectorField& F
) const
{
    if (shape_->centerProc_ == -1 && shape_->singleProc())
    {
        return;
    }

    scalar dT = mesh_.time().deltaT().value();

    const List<labelList>& Is = shape_->Is_;
    const List<labelList>& Os = shape_->Os_;
    const List<scalarList>& ws = shape_->wFromLocal_;
    const scalarList& W = shape_->WFromLocal_;

    forAll(shape_->shellCells_, celli)
    {
        label i = shape_->shellCells_[celli];

        vector Fi = Zero;
        for(label pti = 0; pti < 4; pti++)
        {
            //  Set inner forcing so that u = -u_outer
            Fi +=
                ws[celli][pti]/W[celli]
               *(
                    (
                        v(shape_->baseMesh_[Os[celli][pti]])
                      - U[Os[celli][pti]]
                      - Uold[Is[celli][pti]]
                    )/dT
                  + S[Is[celli][pti]]
                );
        }
        for(label pti = 4; pti < 8; pti++)
        {
            //  Set forcing so that u = u_p
            Fi +=
                ws[celli][pti]/W[celli]
               *(
                    (
                        v(shape_->baseMesh_[Is[celli][pti]])
                      - Uold[Is[celli][pti]]
                    )/dT
                  + S[Is[celli][pti]]
                );
        }
        F[i] += Fi;
    }
}


Foam::vector Foam::particleIBM::v(const vector& pt) const
{
    return (-(pt - shape_->center_)^omega_ ) + v_;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::particleIBM::particleIBM
(
    const polyMesh& mesh,
    const dictionary& dict,
    const label index
)
:
    particle
    (
        (mesh.findCell(dict.lookup("position")) != -1)
      ? particle
        (
            mesh,
            vector(dict.lookup("position")),
            mesh.findCell(dict.lookup("position"))
        )
      : Foam::particle
        (
            mesh,
            mesh.cellCentres()[mesh.findNearestCell(dict.lookup("position"))],
            mesh.findNearestCell(dict.lookup("position"))
        )
    ),
    mesh_(mesh),
    index_(index),
    copy_(0),
    active_(true),
    shape_
    (
        particleShape::New
        (
            mesh,
            dict,
            vector(dict.lookup("position"))
        )
    ),
    v_(dict.template lookupOrDefault<vector>("v", Zero)),
    omega_(dict.template lookupOrDefault<vector>("omega", Zero)),
    rho_(readScalar(dict.lookup("rho"))),
    age_(0.0),
    integratedForce_(Zero),
    integratedForceOld_(Zero),
    integratedTorque_(Zero),
    integratedTorqueOld_(Zero)
{}


Foam::particleIBM::particleIBM
(
    const particleIBM& p,
    const label copy,
    const vector& center,
    const vector& theta,
    const vector& v,
    const vector& omega,
    const scalar& age
)
:
    particle
    (
        (p.mesh_.findCell(center) != -1)
      ? static_cast<const particle&>(p)
      : Foam::particle
        (
            p.mesh_,
            p.mesh_.cellCentres()[p.mesh_.findNearestCell(center)],
            p.mesh_.findNearestCell(center)
        )
    ),
    mesh_(p.mesh_),
    index_(p.index_),
    copy_(copy),
    active_(false),
    shape_(particleShape::New(p.shape_(), center, theta)),
    v_(v),
    omega_(omega),
    rho_(p.rho_),
    age_(age),
    integratedForce_(Zero),
    integratedForceOld_(Zero),
    integratedTorque_(Zero),
    integratedTorqueOld_(Zero),
    collisionForces_(p.collisionForces_.size(), Zero),
    wallForces_(p.wallForces_.size(), Zero)
{}


Foam::particleIBM::particleIBM(const particleIBM& p)
:
    particle
    (
        (p.mesh_.findCell(p.position()) != -1)
      ? static_cast<const particle&>(p)
      : Foam::particle
        (
            p.mesh_,
            p.mesh_.cellCentres()[p.mesh_.findNearestCell(p.position())],
            p.mesh_.findNearestCell(p.position())
        )
    ),
    mesh_(p.mesh_),
    index_(p.index_),
    copy_(p.copy_),
    active_(false),
    shape_(particleShape::New(p.shape_(), p.center(), p.theta())),
    v_(p.v_),
    omega_(p.omega_),
    rho_(p.rho_),
    age_(p.age_),
    integratedForce_(p.integratedForce_),
    integratedForceOld_(p.integratedForceOld_),
    integratedTorque_(p.integratedTorque_),
    integratedTorqueOld_(p.integratedTorqueOld_),
    collisionForces_(p.collisionForces_.size(), Zero),
    wallForces_(p.wallForces_.size(), Zero)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * *//

Foam::particleIBM::~particleIBM()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::particleIBM::solve
(
    const scalar& dt,
    const bool moving,
    const bool rotation
)
{
    if (!onMesh())
    {
        return;
    }

    if (moving)
    {
        vector Fc = sum(collisionForces_);
        vector Fw = sum(wallForces_);
        integratedForce_ += Fc + Fw;

        v_ += dt*integratedForce_/mass();
        vector dx = dt*v_;
        shape_->center_ += dx;
    }

    if (rotation)
    {
        omega_ += dt*integratedTorque_/(mass()*shape_->I());
        shape_->theta() += omega_*dt;
    }
    integratedForce_ = integratedForceOld_;
    integratedTorque_ = integratedTorqueOld_;
    collisionForces_ = Zero;
    wallForces_ = Zero;
}


void Foam::particleIBM::wallHit
(
    const fvMesh& mesh,
    const scalar& dt,
    const label patchi
)
{
    if (!onMesh())
    {
        return;
    }
    forAll(mesh.boundary()[patchi], facei)
    {
        const vector& faceCentre =
            mesh.Cf().boundaryField()[patchi][facei];
        vector norm =
            mesh.Sf().boundaryField()[patchi][facei]
            /mesh.magSf().boundaryField()[patchi][facei];

        vector R = center() - faceCentre;
        scalar magR = mag(R);
        vector normR = R/magR;
        if (magR <= r(faceCentre) && (v_ & norm) > 0)
        {
            vector vNorm = norm*(norm & v_);
            wallForces_[patchi] = -2.0*mass()*(vNorm & normR)*normR/dt;
//                     integratedTorque_ += normalForce_ ^ R;

            // Only calculate one hit per patch
            return;
        }
    }
}

Foam::label Foam::particleIBM::patchHit
(
    const fvMesh& mesh,
    const label patchi,
    vector& R
) const
{
    if (!onMesh())
    {
        return -1;
    }

    forAll(mesh.boundaryMesh()[patchi], facei)
    {
        const vector& faceCentre =
            mesh.Cf().boundaryField()[patchi][facei];
        vector norm =
            mesh.Sf().boundaryField()[patchi][facei]
            /mesh.magSf().boundaryField()[patchi][facei];

        R = center() - faceCentre;
        if (mag(R) <= r(faceCentre) && (v_ & norm) > 0)
        {
            return facei;
        }
    }
    return -1;
}


void Foam::particleIBM::forcing
(
    const surfaceVectorField& Uf,
    const surfaceVectorField& Ufold,
    const surfaceVectorField& Sf,
    volVectorField& F
) const
{
    label N = shape_->N();

    vectorField interpolatedU(N, Zero);
    vectorField interpolatedUold(N, Zero);
    vectorField interpolatedS(N, Zero);

    interpolateFromMesh(Uf,interpolatedU);
    interpolateFromMesh(Ufold,interpolatedUold);
    interpolateFromMesh(Sf,interpolatedS);

    interpolateFToMesh(interpolatedU,interpolatedUold,interpolatedS,F);
}

void Foam::particleIBM::integrateSurfaceStress
(
    const surfaceSymmTensorField& tauf,
    const surfaceScalarField& pf
)
{
    symmTensorField interpolatedTau(shape_->N(), Zero);
    scalarField interpolatedP(shape_->N(), Zero);

    interpolateFromMesh(tauf,interpolatedTau);
    interpolateFromMesh(pf,interpolatedP);
    tmp<vectorField> tmpSf(shape_->Sf());
    const vectorField& Sf = tmpSf();

    integratedForceOld_ = Zero;
    integratedTorqueOld_ = Zero;

    for (label i = 0; i < shape_->nTheta(); i++)
    {
        for (label k = 0; k < shape_->nk(); k++)
        {
            vector R =
                shape_->baseMesh_[shape_->index2(i,k)] - shape_->center_;

            const vector& Sfi = Sf[shape_->index2(i,k)];
            scalar magSf = mag(Sfi);
            vector norm = Sfi/magSf;

            const symmTensor& tau = interpolatedTau[shape_->index(1,i,k)];
            vector F =
                interpolatedP[shape_->index(1,i,k)]
               *shape_->Sf_[shape_->index2(i,k)]
              + (tau & Sfi);
            vector Fn = (F & norm)*norm;
            vector Ft = F - Fn;

            integratedForceOld_ += F;
            integratedTorqueOld_ -= Ft^R;
        }
    }

    integratedForce_ = integratedForceOld_;
    integratedTorque_ = integratedTorqueOld_;
}


Foam::scalar Foam::particleIBM::Cd
(
    const scalar& rhoRef,
    const vector& UInf
) const
{
    vector normUInf(UInf/mag(UInf));

    dimensionedScalar Cd
    (
        "Cd",
        dimless,
        (normUInf & integratedForce_)*2.0
        /(rhoRef*magSqr(UInf)*shape_->A(shape_->center_ - normUInf))
    );

    return Cd.value();
}

void Foam::particleIBM::update()
{
    shape_->moveMesh(shape_->center_);
    shape_->updateCellLists();

    if (shape_->centerOnMesh())
    {
        active_ = 1;
    }
    else
    {
        active_ = 0;
    }

    if (active_)
    {
        track(shape_->center() - position(), 1.0);
    }
}