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

#include "particle.H"
#include "particleShape.H"
#include "particleList.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::IBM::particle::interpolateFToMesh
(
    const vectorList& U,
    const vectorList& Uold,
    const vectorList& S,
    volVectorField& F
)
{
    if (origProc_ == -1 && neiProcs_.size() == 0)
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


void Foam::IBM::particle::setProcs()
{
    if (mesh_.findCell(shape_->CoM()) != -1)
    {
        origProc_ = Pstream::myProcNo();
    }

    neiProcs_.clear();

    forAll(shape_->shellCells(), celli)
    {
        if
        (
            shape_->shellCells()[celli] != -1
         && Pstream::myProcNo() != origProc_
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

Foam::IBM::particle::particle
(
    const particleList<particle>& pList,
    const dictionary& dict,
    const volVectorField& Uold,
    const volVectorField& S
)
:
    pList_(pList),
    mesh_(pList_.mesh()),
    shape_(particleShape::New(mesh_,dict)),
    U_(mesh_.lookupObject<volVectorField>("U")),
    U0_(Uold),
//     rho_(mesh_.lookupObject<volScalarField>("thermo:rho")),
    p_(mesh_.lookupObject<volScalarField>("p")),
    S_(S),
//     UInterp_
//     (
//         interpolation<vector>::New
//         (
//             pList.ibmDict().subDict("interpolationSchemes"),
//             U_
//         )
//     ),
//     U0Interp_
//     (
//         interpolation<vector>::New
//         (
//             pList.ibmDict().subDict("interpolationSchemes"),
//             U0_
//         )
//     ),
//     pInterp_
//     (
//         interpolation<scalar>::New
//         (
//             pList.ibmDict().subDict("interpolationSchemes"),
//             p_
//         )
//     ),
//     SInterp_
//     (
//         interpolation<vector>::New
//         (
//             pList.ibmDict().subDict("interpolationSchemes"),
//             S_
//         )
//     ),
    v_(dict.template lookupOrDefault<vector>("v", Zero)),
    omega_(dict.template lookupOrDefault<vector>("omega", Zero)),
    rho_(readScalar(dict.lookup("rho"))),
    integratedForce_(Zero),
    origProc_(-1)
{
    setProcs();
}


Foam::IBM::particle::~particle()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

bool Foam::IBM::particle::onMesh() const
{
    if (shape_->shellCells().size() == 0)
    {
        return false;
    }
    return true;
}

void Foam::IBM::particle::forcing
(
    const surfaceVectorField& Uf,
    const surfaceVectorField& Ufold,
    const surfaceVectorField& Sf,
    volVectorField& F
)
{
    label N = shape_->N();

    vectorList interpolatedU(N, Zero);
    vectorList interpolatedUold(N, Zero);
    vectorList interpolatedS(N, Zero);

    interpolateFromMesh(Uf,interpolatedU);
    interpolateFromMesh(Ufold,interpolatedUold);
    interpolateFromMesh(Sf,interpolatedS);

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

void Foam::IBM::particle::integrateSurfaceStress
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

Foam::scalar Foam::IBM::particle::Cd() const
{
    dimensionedScalar Cd
        (
            "Cd",
            dimless,
            (pList_.normUInf() & integratedForce_)*2.0
           /(pList_.rhoRef()*magSqr(pList_.UInf())*shape_->A())
        );

    combineReduce
    (
        Cd,
        MaxEqOp<dimensionedScalar>()
    );
    return Cd.value();
}

void Foam::IBM::particle::update()
{
    shape_->moveMesh();
    shape_->updateCellLists();
    setProcs();
    shape_->updateCellLists();
}