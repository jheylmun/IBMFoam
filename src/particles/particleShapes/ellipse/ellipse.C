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

#include "ellipse.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::particleShapes::ellipse::calcSf()
{
    scalar Pi = constant::mathematical::pi;
    scalar dTheta = 2.0*Pi/scalar(nTheta_);

    for (label i = 0; i < nTheta_; i++)
    {
        Sf_[index2(i,0)] = -centeredMesh_[index(1,i,0)]*dTheta*l_;
    }
}

Foam::tmp<Foam::vectorField> Foam::particleShapes::ellipse::rotate()
{
    return centeredMesh_ & rotationMatrix_;
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleShapes
{
    defineTypeNameAndDebug(ellipse, 0);

    addToRunTimeSelectionTable
    (
        particleShape,
        ellipse,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleShapes::ellipse::ellipse
(
    const polyMesh& mesh,
    const dictionary& dict,
    const vector& center
)
:
    particleShape(mesh, dict, center),
    a_(readScalar(dict.lookup("a"))),
    b_(readScalar(dict.lookup("b")))
{
    if (mesh.nGeometricD() != 2)
    {
        FatalErrorInFunction
            << "ellipse particle shape with 1 point in axis " << nl
            << "only valid for 2D meshes"
            << exit(FatalError);
    }

    this->momentOfInertia_ = a_*b_/4.0;
    nk_ = 1;

    boundBox bb(mesh.points());

    l_ = mag(bb.max().z() - bb.min().z());

    discretize();
    updateCellLists();
    calcSf();
}

Foam::particleShapes::ellipse::ellipse
(
    const particleShape& shape,
    const vector& center,
    const vector& theta
)
:
    particleShape(shape, center, theta),
    a_(refCast<const ellipse>(shape).a_),
    b_(refCast<const ellipse>(shape).b_),
    l_(refCast<const ellipse>(shape).l_)
{
    discretize();
    updateCellLists();
    calcSf();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleShapes::ellipse::~ellipse()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::particleShapes::ellipse::discretize()
{
    scalar Pi = constant::mathematical::pi;
    vector delta(-1*delta_, 0, delta_);
    scalar dTheta(2.0*Pi/(nTheta_));

    for (label i = 0; i < nRadial_; i++)
    {
        for (label j = 0; j < nTheta_; j++)
        {
            scalar theta = scalar(j)*dTheta;
            scalar R =
                sqrt(1.0/(sqr(cos(theta)/a_) + sqr(sin(theta)/b_)))
              + delta[i];
            label celli = index(i,j,0);

            centeredMesh_[celli] =
                vector
                (
                    R*Foam::cos(theta),
                    R*Foam::sin(theta),
                    scalar(0)
                );
        }
    }
    this->moveMesh(this->center());
}


void Foam::particleShapes::ellipse::updateCellLists()
{
    scalar twoPi = constant::mathematical::twoPi;
    if (theta_.z() < 0) theta_.z() += twoPi;
    if (theta_.z() > twoPi) theta_.z() -= twoPi;

    shellCells_ = labelList(mesh_.nCells(), -1);
    neighbourPoints_ = List<labelVector>(mesh_.nCells(), Zero);

    label i = 0;
    forAll(mesh_.cellCentres(), celli)
    {
        const vector& CC = mesh_.cellCentres()[celli];
        vector diff = CC - center_;

        scalar radius = mag(diff);
        scalar R = r(CC);
        scalar innerR = R - delta_;

        if (radius >= innerR && radius <= R)
        {
            shellCells_[i] = celli;

            scalar theta = Foam::atan2(diff.y(),diff.x()) - theta_.z();
            if (theta < 0) theta += twoPi;

            neighbourPoints_[i] =
            (
                labelVector
                (
                    0,
                    label(theta*(nTheta_ - 1)/twoPi),
                    0
                )
            );
            i++;
        }
    }
    shellCells_.resize(i);
    neighbourPoints_.resize(i);

    this->setNeighbours();
    this->setWeights();
    setProcs();
}

Foam::scalar Foam::particleShapes::ellipse::d() const
{
    return a_*b_*Foam::constant::mathematical::twoPi*sqrt((sqr(a_) + sqr(b_))/2.0);
}

Foam::scalar Foam::particleShapes::ellipse::r(const vector& pt) const
{
    vector diff = pt - center_;
    scalar theta = Foam::atan2(diff.y(), diff.x()) - theta_.z();
    return sqrt(1.0/(sqr(cos(theta)/a_) + sqr(sin(theta)/b_)));
}

Foam::scalar Foam::particleShapes::ellipse::A(const vector& pt) const
{
    return a_*b_*Foam::constant::mathematical::pi;
}

Foam::scalar Foam::particleShapes::ellipse::V() const
{
    return
        Foam::constant::mathematical::pi*a_*b_*l_/4.0;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::particleShapes::ellipse::write
(
    Ostream& os,
    const particleShape& shape
) const
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particleShapes::ellipse&>(shape).d();
    }
    else
    {
        os  << static_cast<const particleShapes::ellipse&>(shape);
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const particleShape&)"
    );

    return os;
}

// ************************************************************************* //