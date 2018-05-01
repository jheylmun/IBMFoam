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
namespace particleShapes
{
    defineTypeNameAndDebug(cylinder, 0);

    addToRunTimeSelectionTable
    (
        particleShape,
        cylinder,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        particleShape,
        cylinder,
        copy
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleShapes::cylinder::cylinder
(
    const polyMesh& mesh,
    const dictionary& dict,
    const vector& center
)
:
    particleShape(mesh, dict, center),
    d_(readScalar(dict.lookup("d")))
{
    if (mesh.nGeometricD() != 2)
    {
        FatalErrorInFunction
            << "Cylinder particle shape with 1 point in axis " << nl
            << "only valid for 2D meshes"
            << exit(FatalError);
    }
    momentOfInertia_ = Foam::sqr(d_/8.0);
    nk_ = 1;
    l_ = mag(max(mesh.points()).z() - min(mesh.points()).z());

    discretize();
    updateCellLists();
    calcSf();

}

Foam::particleShapes::cylinder::cylinder
(
    const particleShape& shape,
    const vector& center,
    const vector& theta
)
:
    particleShape(shape, center, theta),
    d_(refCast<const cylinder>(shape).d_),
    l_(refCast<const cylinder>(shape).l_)
{
    this->moveMesh(center);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleShapes::cylinder::~cylinder()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::particleShapes::cylinder::calcSf()
{
    scalar r = d_/2.0;
    scalar Pi = constant::mathematical::pi;
    scalar dTheta = 2.0*Pi/scalar(nTheta_);

    for (label i = 0; i < nTheta_; i++)
    {
        scalar theta = dTheta*scalar(i);
        Sf_[index2(i,0)] =
           -vector
            (
                cos(theta),
                sin(theta),
                scalar(0)
            )*r*dTheta*l_;
    }
}

void Foam::particleShapes::cylinder::discretize()
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

            label celli = index(i,j,0);

            centeredMesh_[celli] =
                vector
                (
                    r*Foam::cos(theta),
                    r*Foam::sin(theta),
                    scalar(0)
                );
        }
    }
    this->moveMesh(this->center());
}


void Foam::particleShapes::cylinder::updateCellLists()
{
    scalar R = d_/2.0;
    scalar innerR = R - delta_;
    scalar twoPi = constant::mathematical::twoPi;

    shellCells_ = labelList(mesh_.nCells(), -1);
    neighbourPoints_ = List<labelVector>(mesh_.nCells(), Zero);

    label i = 0;
    forAll(mesh_.cellCentres(), celli)
    {
        vector diff = mesh_.cellCentres()[celli] - center_;

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
}

Foam::scalar Foam::particleShapes::cylinder::d() const
{
    return d_;
}


Foam::vector Foam::particleShapes::cylinder::D() const
{
    return vector(d_, d_, d_);
}


Foam::scalar Foam::particleShapes::cylinder::A(const vector& pt) const
{
    return d_*l_;
}

Foam::scalar Foam::particleShapes::cylinder::V() const
{
    return
        Foam::constant::mathematical::pi
       *sqr(d_/2.0)*l_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::particleShapes::cylinder::write
(
    Ostream& os,
    const particleShape& shape
) const
{
    const cylinder& c = static_cast<const particleShapes::cylinder&>(shape);
    os  << "d" << token::TAB << c.d_ << token::END_STATEMENT << endl;

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const particleShapes::ellipse&)"
    );

    return os;
}

// ************************************************************************* //