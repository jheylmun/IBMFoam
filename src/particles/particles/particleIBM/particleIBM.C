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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(particleIBM, 0);
    defineTemplateTypeNameAndDebug(Cloud<particleIBM>, 0);
}

Foam::label Foam::particleIBM::centerOfMesh(const polyMesh& mesh)
{
    vector center(0.5*(max(mesh.points()) + min(mesh.points())));
    return mesh.findNearestCell(center);
}

bool Foam::particleIBM::collision
(
    const particleIBM& p1,
    const particleIBM& p2,
    vector& collisionPt
)
{
    scalar pi = Foam::constant::mathematical::pi;

    vector center1 = p1.center();
    vector center2 = p2.center();
    vector diff = center2 - center1;

    vector d1 = p1.D();
    vector d2 = p2.D();
    if (mag(diff) > 0.5*(cmptMax(d1) + cmptMax(d2)))
    {
        return false;
    }
    if
    (
        (d1.x() + d1.y() + d1.z())/3.0 == cmptMax(d1)
     && (d2.x() + d2.y() + d2.z())/3.0 == cmptMax(d2)
    )
    {
        if (0.5*(cmptMax(d1) + cmptMax(d2)) > mag(diff))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    vector pt1 = center1;
    vector pt2 = center2;
    tensor R1 = p1.rotationMatrix();
    tensor R2 = p2.rotationMatrix();

    scalar gamma1 = atan2(diff.y(), diff.x()) - p1.theta().z();
    scalar gamma2 = atan2(diff.y(), diff.x()) - p2.theta().z();
    if (center2.x() < center1.x())
    {
        gamma2 -= pi;
    }
    else
    {
        gamma1 -= pi;
    }

    scalar r1 = p1.r(pt2);
    scalar r2 = p2.r(pt1);
    scalar x1 = r1*cos(gamma1);
    scalar x2 = r2*cos(gamma2);
    if (tan(p1.theta().z())*tan(p2.theta().z()) < 0)
    {
        x1 *= -1.0;
        x2 *= -1.0;
    }

    bool converged = false;
    bool collision = false;
    label itt = 0;

    while (!converged)
    {
        itt++;
        pt1 = (R1 & vector(x1, 0, 0)) + center1;
        pt2 = (R2 & vector(x2, 0, 0)) + center2;

        vector oldDiff = diff;
        diff = pt2 - pt1;

        scalar alpha1 = atan2(diff.y(), diff.x()) - p1.theta().z();
        scalar alpha2 = atan2(diff.y(), diff.x()) - p2.theta().z();
        if (center2.x() < center1.x())
        {
            alpha2 -= pi;
        }
        else
        {
            alpha1 -= pi;
        }

        vector rp1 = pt1;
        vector rp1Local = rp1;
        {
            scalar m = tan(alpha1);
            scalar I = -m*x1;
            scalar a = sqr(m) + sqr(d1.y()/d1.x());
            scalar b = 2.0*m*I;
            scalar c = sqr(I) - sqr(d1.y());

            scalar posRoot = (-b + sqrt(sqr(b) - 4.0*a*c))/(2.0*a);
            scalar negRoot = (-b - sqrt(sqr(b) - 4.0*a*c))/(2.0*a);

            vector rpPos = (R1 & vector(posRoot, I + posRoot*m, 0)) + center1;
            vector rpNeg = (R1 & vector(negRoot, I + negRoot*m, 0)) + center1;

            if (mag(rpPos - center2) < mag(rpNeg - center2))
            {
                rp1 = rpPos;
                rp1Local = vector(posRoot, I + posRoot*m, 0);
            }
            else
            {
                rp1 = rpNeg;
                rp1Local = vector(negRoot, I + negRoot*m, 0);
            }
            gamma1 = atan2(rp1Local.y(), rp1Local.x());
        }

        vector rp2 = pt2;
        vector rp2Local = rp2;
        {
            scalar m = tan(alpha2);
            scalar I = -m*x2;
            scalar a = sqr(m) + sqr(d2.y()/d2.x());
            scalar b = 2.0*m*I;
            scalar c = sqr(I) - sqr(d2.y());

            scalar posRoot = (-b + sqrt(sqr(b) - 4.0*a*c))/(2.0*a);
            scalar negRoot = (-b - sqrt(sqr(b) - 4.0*a*c))/(2.0*a);

            vector rpPos = (R2 & vector(posRoot, I + posRoot*m, 0)) + center2;
            vector rpNeg = (R2 & vector(negRoot, I + negRoot*m, 0)) + center2;

            if (mag(rpPos - center1) < mag(rpNeg - center1))
            {
                rp2 = rpPos;
                rp2Local = vector(posRoot, I + posRoot*m, 0);
            }
            else
            {
                rp2 = rpNeg;
                rp2Local = vector(negRoot, I + negRoot*m, 0);
            }
            gamma2 = atan2(rp2Local.y(), rp2Local.x());
        }

        r1 = mag(rp1Local);
        r2 = mag(rp2Local);
        x1 = r1*cos(gamma1);
        x2 = r2*cos(gamma2);
        if (tan(p1.theta().z())*tan(p2.theta().z()) < 0)
        {
            x1 *= -1.0;
            x2 *= -1.0;
        }

        if ((diff & (rp2 - rp1)) < 0)
        {
            converged = true;
            collision = true;
            collisionPt = 0.5*(rp1 + rp2);
        }

        if (mag(oldDiff) < mag(diff) || itt > 10)
        {
            converged = true;
        }
    }
    return collision;
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
    if (!onMesh())
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
                        v(shape_->mesh_[Os[celli][pti]])
                      - U[Os[celli][pti]]
                      - Uold[Is[celli][pti]]
                    )/dT
                  + S[Os[celli][pti]]
                );
        }
        for(label pti = 4; pti < 8; pti++)
        {
            //  Set forcing so that u = u_p
            Fi +=
                ws[celli][pti]/W[celli]
               *(
                    (
                        v(shape_->mesh_[Is[celli][pti]])
                      - Uold[Is[celli][pti]]
                    )/dT
                  + S[Is[celli][pti]]
                );
        }
        F[i] += Fi;
    }
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
        mesh,
        mesh.cellCentres()[mesh.findNearestCell(dict.lookup("position"))],
        mesh.findNearestCell(dict.lookup("position"))
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
    integratedTorque_(Zero)
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
        p.mesh_,
        p.mesh_.cellCentres()[centerOfMesh(p.mesh_)],
        centerOfMesh(p.mesh_)
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
    integratedForce_(p.integratedForce_),
    integratedTorque_(p.integratedTorque_),
    collisionForces_(p.collisionForces_.size(), Zero),
    collisionTorques_(p.collisionTorques_.size(), Zero),
    wallForces_(p.wallForces_.size(), Zero),
    wallTorques_(p.wallTorques_.size(), Zero)
{}


Foam::particleIBM::particleIBM(const particleIBM& p)
:
    particle
    (
        p.mesh_,
        p.mesh_.cellCentres()[centerOfMesh(p.mesh_)],
        centerOfMesh(p.mesh_)
    ),
    mesh_(p.mesh_),
    index_(p.index_),
    copy_(p.copy_),
    active_(false),
    shape_(),//particleShape::New(p.shape_(), p.center(), p.theta())),
    v_(p.v_),
    omega_(p.omega_),
    rho_(p.rho_),
    age_(p.age_),
    integratedForce_(p.integratedForce_),
    integratedTorque_(p.integratedTorque_),
    collisionForces_(p.collisionForces_.size(), Zero),
    collisionTorques_(p.collisionTorques_.size(), Zero),
    wallForces_(p.wallForces_.size(), Zero),
    wallTorques_(p.wallTorques_.size(), Zero)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * *//

Foam::particleIBM::~particleIBM()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

Foam::vector Foam::particleIBM::v(const vector& pt) const
{
    return (-(pt - shape_->center_)^omega_ ) + v_;
}


void Foam::particleIBM::solve
(
    const scalar& dt,
    const vector& g,
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

    if (rotation)
    {
        vector T =
            integratedTorque_
          + sum(collisionTorques_)
          + sum(wallTorques_);

        omega_ += dt*T/(mass()*shape_->I());
        shape_->theta() += omega_*dt;

        scalar twoPi = constant::mathematical::twoPi;
        if (shape_->theta().z() < 0) shape_->theta().z() += twoPi;
        if (shape_->theta().z() >= twoPi) shape_->theta().z() -= twoPi;
    }
    collisionForces_ = Zero;
    collisionTorques_ = Zero;
    wallForces_ = Zero;
    wallTorques_ = Zero;
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

        vector rv = v(faceCentre);

        vector R = center() - faceCentre;
        if (mag(R) <= r(faceCentre) && (rv & norm) > 0)
        {
            vector vNorm = norm*(norm & rv);

            vector rHat = R/mag(R);
            vector F = -2.0*mass()*vNorm/dt;
            wallForces_[patchi] = (F & rHat)*rHat;
            wallTorques_[patchi] = (F - wallForces_[patchi])^R;

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

    integratedForce_ = Zero;
    integratedTorque_ = Zero;

    for (label i = 0; i < shape_->nTheta(); i++)
    {
        for (label k = 0; k < shape_->nk(); k++)
        {
            vector R =
                shape_->mesh_[shape_->index2(i,k)] - shape_->center_;

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

            integratedForce_ += F;
            integratedTorque_ -= Ft^R;
        }
    }
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