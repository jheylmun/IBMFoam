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
    if (centerProc_ == -1 && neiProcs_.size() == 0)
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
                        this->v() - U[Os[celli][pti]]
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
                        this->v()
                      - Uold[Is[celli][pti]]
                    )/dT
                  + S[Is[celli][pti]]
                );
        }
        F[i] += Fi;
    }
}


void Foam::particleIBM::setProcs()
{
    if (mesh_.findCell(shape_->center()) != -1)
    {
        centerProc_ = Pstream::myProcNo();
    }

    neiProcs_.clear();

    forAll(shape_->shellCells(), celli)
    {
        if
        (
            shape_->shellCells()[celli] != -1
         && Pstream::myProcNo() != centerProc_
        )
        {
            bool set = false;
            forAll(neiProcs_, proci)
            {
                if (neiProcs_[proci] == Pstream::myProcNo())
                {
                    set = true;
                    break;
                }
            }
            if (!set)
            {
                neiProcs_.append(Pstream::myProcNo());
            }
        }
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::particleIBM::particleIBM
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    particle
    (
        mesh,
        vector(dict.lookup("position")),
        mesh.findCell(dict.lookup("position"))
    ),
    mesh_(mesh),
    dict_(dict),
    active_(true),
    shape_(particleShape::New(mesh_, dict, this->position())),
    v_(dict.template lookupOrDefault<vector>("v", Zero)),
    omega_(dict.template lookupOrDefault<vector>("omega", Zero)),
    rho_(readScalar(dict.lookup("rho"))),
    age_(0.0),
    integratedForce_(Zero),
    centerProc_(-1)
{
    setProcs();
}


Foam::particleIBM::particleIBM
(
    const polyMesh& mesh,
    const dictionary& dict,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    particle(mesh, coordinates, celli, tetFacei, tetPti),
    mesh_(mesh),
    dict_(dict),
    active_(true),
    shape_(particleShape::New(mesh_, dict, position())),
    v_(dict.template lookupOrDefault<vector>("v", Zero)),
    omega_(dict.template lookupOrDefault<vector>("omega", Zero)),
    rho_(readScalar(dict.lookup("rho"))),
    age_(0.0),
    integratedForce_(Zero),
    centerProc_(-1)
{
    setProcs();
}


Foam::particleIBM::particleIBM
(
    const polyMesh& mesh,
    const dictionary& dict,
    const vector& position,
    const label celli
)
:
    particle(mesh, position, celli),
    mesh_(mesh),
    dict_(dict),
    active_(true),
    shape_(particleShape::New(mesh_, dict_, this->position())),
    v_(dict.template lookupOrDefault<vector>("v", Zero)),
    omega_(dict.template lookupOrDefault<vector>("omega", Zero)),
    rho_(readScalar(dict.lookup("rho"))),
    age_(0.0),
    integratedForce_(Zero),
    centerProc_(-1)
{
    setProcs();
}


Foam::particleIBM::particleIBM
(
    const particleIBM& p,
    const polyMesh& mesh
)
:
    particle(p),
    mesh_(mesh),
    dict_(p.dict_),
    active_(p.active_),
    shape_(particleShape::New(mesh_, dict_, position())),
    v_(p.v_),
    omega_(p.omega_),
    rho_(p.rho_),
    age_(p.age_),
    integratedForce_(p.integratedForce_),
    centerProc_(p.centerProc_)
{
    setProcs();
}


Foam::particleIBM::particleIBM
(
    const particleIBM& p
)
:
    particle(p),
    mesh_(p.mesh_),
    dict_(p.dict_),
    active_(p.active_),
    shape_(particleShape::New(mesh_, dict_, particle::position())),
    v_(p.v_),
    omega_(p.omega_),
    rho_(p.rho_),
    age_(p.age_),
    integratedForce_(p.integratedForce_),
    centerProc_(p.centerProc_)
{
    setProcs();
}


Foam::particleIBM::~particleIBM()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::particleIBM::solve
(
    const scalar& dt
)
{
    vector disp = dt*v_;

    shape_->center_ += disp;
    v_ += dt*integratedForce_/mass();
    update();

    this->track(disp, 1.0);
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

    vectorList interpolatedU(N, Zero);
    vectorList interpolatedUold(N, Zero);
    vectorList interpolatedS(N, Zero);

    interpolateFromMesh<vector>(Uf,interpolatedU);
    interpolateFromMesh<vector>(Ufold,interpolatedUold);
    interpolateFromMesh<vector>(Sf,interpolatedS);

    if (Pstream::parRun())
    {
        combineReduce
        (
            interpolatedU,
            ListPlusEqOp<vectorList>()
        );
        combineReduce
        (
            interpolatedUold,
            ListPlusEqOp<vectorList>()
        );
        combineReduce
        (
            interpolatedS,
            ListPlusEqOp<vectorList>()
        );
    }

    interpolateFToMesh(interpolatedU,interpolatedUold,interpolatedS,F);
}

void Foam::particleIBM::integrateSurfaceStress
(
    const surfaceSymmTensorField& tauf,
    const surfaceScalarField& pf
)
{
    symmTensorList interpolatedTau(shape_->N(), Zero);
    scalarList interpolatedP(shape_->N(), Zero);

    interpolateFromMesh(tauf,interpolatedTau);
    interpolateFromMesh(pf,interpolatedP);
    if (Pstream::parRun())
    {
        combineReduce
        (
            interpolatedTau,
            ListPlusEqOp<symmTensorList>()
        );
        combineReduce
        (
            interpolatedP,
            ListPlusEqOp<scalarList>()
        );
    }

    integratedForce_ = Zero;

    for (label i = 0; i < shape_->nTheta(); i++)
    {
        for (label k = 0; k < shape_->nk(); k++)
        {
            integratedForce_ +=
                interpolatedP[shape_->index(1,i,k)]
               *shape_->Sf_[shape_->index2(i,k)]
              + (
                    interpolatedTau[shape_->index(1,i,k)]
                  & shape_->Sf_[shape_->index2(i,k)]
                );
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
           /(rhoRef*magSqr(UInf)*shape_->A())
        );

    combineReduce
    (
        Cd,
        MaxEqOp<dimensionedScalar>()
    );
    return Cd.value();
}

void Foam::particleIBM::update()
{
    shape_->moveMesh(shape_->center_);
    shape_->updateCellLists();
    setProcs();
    shape_->updateCellLists();

    if (centerProc_ != -1)
    {
        position() = shape_->center_;
    }
}