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

#include "cylinder.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace IBM
{
namespace particleShapes
{
    defineTypeNameAndDebug(cylinder, 0);

    addToRunTimeSelectionTable
    (
        particleShape,
        cylinder,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IBM::particleShapes::cylinder::cylinder
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    particleShape(mesh,dict),
    d_(readScalar(dict.lookup("d"))),
    p1_(dict.lookup("p1")),
    p2_(dict.lookup("p2")),
    centerPoint_((p1_ + p2_)/2.0),
    axis_(nk_, Zero),
    centeredAxis_(nk_, Zero)
{
    if (!dict_.found("nk_"))
    {
        nk_ = 1;
    }

    if (nk_ == 1)
    {
        if (mesh.nGeometricD() != 2)
        {
            FatalErrorInFunction
                << "Cylinder particle shape with 1 point in axis " << nl
                << "only valid for 2D meshes"
                << exit(FatalError);
        }

        axis_[0] = centerPoint_;
        centeredAxis_[0] - axis_[0] - CoM();
    }
    else
    {
        forAll(axis_, k)
        {
            axis_[k] = p1_ + scalar(k)/(nk_ - 1)*(p2_ - p1_);
            centeredAxis_[k] = axis_[k] - CoM();
        }
    }

    discretize();
    updateCellLists();
    calcSf();
}

void Foam::IBM::particleShapes::cylinder::calcSf()
{
    scalar r = d_/2.0;
    scalar Pi = constant::mathematical::pi;
    scalar dTheta = 2.0*Pi/scalar(nTheta_);
    scalar dk = mag(p1_ - p2_)/(max(1,nk_ - 1));

    for (label i = 0; i < nTheta_; i++)
    {
        scalar theta = dTheta*scalar(i);
        for (label k = 0; k < nk_; k++)
        {
            Sf_[index2(i,k)] =
               -vector
                (
                    cos(theta),
                    sin(theta),
                    scalar(0)
                )*r*dTheta*dk;
            if (nk_ != 1 && (k == 0 || k == nk_ - 1))
            {
                Sf_[index2(i,k)] /= 2.0;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IBM::particleShapes::cylinder::~cylinder()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IBM::particleShapes::cylinder::discretize()
{
    scalar Pi = constant::mathematical::pi;
    scalar R = d_/2.0;
    vector delta(-1*delta_, 0, delta_);
    scalar dTheta(2.0*Pi/(nTheta_));

    for (label i = 0; i < nRadial_; i++)
    {
        scalar r = R + delta[i];
        for (label j = 0; j < nTheta_; j++)
        {
            scalar theta = scalar(j)*dTheta;

            for (label k = 0; k < nk_; k++)
            {
                label celli = index(i,j,k);

                centeredMesh_[celli] =
                    vector
                    (
                        r*Foam::cos(theta),
                        r*Foam::sin(theta),
                        scalar(0)
                    );
            }
        }
    }
    this->moveMesh();
}


void Foam::IBM::particleShapes::cylinder::updateCellLists()
{
    scalar R = d_/2.0;
    scalar innerR = R - delta_;
    scalar twoPi = constant::mathematical::twoPi;

    shellCells_ = labelList(mesh_.nCells(), -1);
    neighbourPoints_ = List<labelVector>(mesh_.nCells(), Zero);

    label i = 0;
    forAll(mesh_.cellCentres(), celli)
    {
        for (label k = 0; k < nk_; k++)
        {
            vector diff = mesh_.cellCentres()[celli] - axis_[k];

            scalar r = mag(diff);

            if (r >= innerR && r <= R)
            {
                shellCells_[i] = celli;

                scalar theta = Foam::atan2(diff.y(),diff.x());

                if (theta < 0) theta += twoPi;

                neighbourPoints_[i] =
                (
                    labelVector
                    (
                        0,
                        label(theta*(nTheta_ - 1)/twoPi),
                        k
                    )
                );
                i++;
            }
        }
    }
    shellCells_.resize(i);
    neighbourPoints_.resize(i);

    this->setNeighbours();
    this->setWeights();
}

Foam::scalar Foam::IBM::particleShapes::cylinder::D() const
{
    return d_;
}

Foam::scalar Foam::IBM::particleShapes::cylinder::A() const
{
    return d_*mag(p2_ - p1_);
}

Foam::scalar Foam::IBM::particleShapes::cylinder::V() const
{
    return
        Foam::constant::mathematical::pi
       *sqr(d_/2.0)*mag(p1_.z() - p2_.z());
}

const Foam::vector& Foam::IBM::particleShapes::cylinder::CoM() const
{
    return centerPoint_;
}

Foam::vector& Foam::IBM::particleShapes::cylinder::CoM()
{
    return centerPoint_;
}

void Foam::IBM::particleShapes::cylinder::moveMesh()
{
    forAll(baseMesh_, celli)
    {
        baseMesh_[celli] = centeredMesh_[celli] + CoM();
    }
    forAll(axis_, k)
    {
        axis_[k] = centeredAxis_[k] + CoM();
    }
}