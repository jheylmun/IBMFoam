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

#include "sphere.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleShapes
{
    defineTypeNameAndDebug(sphere, 0);

    addToRunTimeSelectionTable
    (
        particleShape,
        sphere,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        particleShape,
        sphere,
        copy
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleShapes::sphere::sphere
(
    const polyMesh& mesh,
    const dictionary& dict,
    const vector& center,
    const bool buildMesh
)
:
    particleShape(mesh, dict, center, buildMesh),
    d_(readScalar(dict.lookup("d")))
{
    this->momentOfInertia_ = 2.0/5.0*sqr(d_/2.0);

    if (!dict.found("nk_"))
    {
        nk_ = nTheta_/2;
        N_ = nRadial_*nTheta_*nk_;
        centeredMesh_ = vectorField(N_, Zero);
        mesh_ = vectorField(N_, Zero);
        Sf_ = vectorField(nTheta_*nk_, Zero);
        neighbourPoints_ = List<labelVector>(N_);
        wToLocal_ = List<scalarList>(N_);
        WToLocal_ = scalarList(N_);
        facesToLocal_ = List<labelList>(N_);
    }

    if (mesh.nGeometricD() != 3)
    {
        FatalErrorInFunction
            << "Sphere particle shape on valid for 3D meshes"
            << exit(FatalError);
    }

    if (this->buildMesh_)
    {
        discretize();
        updateCellLists();
        calcSf();
    }
}


Foam::particleShapes::sphere::sphere
(
    const particleShape& shape,
    const vector& center,
    const vector& theta,
    const bool buildMesh
)
:
    particleShape(shape, center, theta, buildMesh),
    d_(refCast<const sphere>(shape).d_)
{
    this->moveMesh(center);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleShapes::sphere::~sphere()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::particleShapes::sphere::calcSf()
{
    scalar r = d_/2.0;
    scalar Pi = constant::mathematical::pi;
    scalar dTheta = 2.0*Pi/scalar(nTheta_);
    scalar dPhi = Pi/scalar(nk_ - 1);

    for (label i = 0; i < nTheta_; i++)
    {
        scalar theta = dTheta*scalar(i);

        for (label k = 0; k < nk_; k++)
        {
            scalar phi = dPhi*scalar(k);

            Sf_[index2(i,k)] =
               -vector
                (
                    Foam::cos(theta)*Foam::sin(phi),
                    Foam::sin(theta)*Foam::sin(phi),
                    Foam::cos(phi)
                )*Foam::sqr(r)*dTheta;

            if (k == 0)
            {
                Sf_[index2(i,k)] *= (cos(phi) - cos(phi + dPhi/2.0));
            }
            else if (k == nk_ - 1)
            {
                Sf_[index2(i,k)] *= (cos(phi - dPhi/2.0) - cos(phi));
            }
            else
            {
                Sf_[index2(i,k)] *=
                    (
                        cos(phi - dPhi/2.0)
                      - cos(phi + dPhi/2.0)
                    );
            }
        }
    }
}


void Foam::particleShapes::sphere::discretize()
{
    scalar Pi = constant::mathematical::pi;
    scalar R = d_/2.0;
    vector delta(-1*delta_, 0, delta_);
    scalar dTheta(2.0*Pi/(nTheta_));
    scalar dPhi = Pi/(nk_ - 1);

    for (label i = 0; i < nRadial_; i++)
    {
        scalar r = R + delta[i];
        for (label j = 0; j < nTheta_; j++)
        {
            scalar theta = scalar(j)*dTheta;

            for (label k = 0; k < nk_; k++)
            {
                label celli = index(i,j,k);

                scalar phi = dPhi*scalar(k);

                centeredMesh_[celli] =
                    vector
                    (
                        r*Foam::cos(theta)*Foam::sin(phi),
                        r*Foam::sin(theta)*Foam::sin(phi),
                        r*Foam::cos(phi)
                    );
            }
        }
    }
    this->moveMesh(this->center_);
}


void Foam::particleShapes::sphere::updateCellLists()
{
    scalar R = d_/2.0;
    scalar innerR = R - delta_;
    scalar Pi = constant::mathematical::pi;

    shellCells_ = labelList(pMesh_.nCells(), -1);
    neighbourPoints_ = List<labelVector>(pMesh_.nCells(), Zero);

    label i = 0;
    forAll(pMesh_.cellCentres(), celli)
    {
        vector diff =
            pMesh_.cellCentres()[celli]
          - this->center_;

        scalar r = mag(diff);

        if (r >= innerR && r <= R)
        {
            shellCells_[i] = celli;

            scalar theta = Foam::atan2(diff.y(),diff.x());
            scalar phi = Foam::acos(diff.z()/r);

            if (theta < 0) theta += 2.0*Pi;
            if (phi < 0) phi += 2.0*Pi;

            neighbourPoints_[i] =
            (
                labelVector
                (
                    0,
                    label(theta*nTheta_/(2.0*Pi)),
                    label(phi*nk_/Pi)
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


Foam::scalar Foam::particleShapes::sphere::d() const
{
    return d_;
}

Foam::vector Foam::particleShapes::sphere::D() const
{
    return vector(d_, d_, d_);
}

Foam::scalar Foam::particleShapes::sphere::A(const vector& pt) const
{
    return Foam::constant::mathematical::pi*sqr(d_/2.0);
}


Foam::scalar Foam::particleShapes::sphere::V() const
{
    return 1.0/6.0*Foam::constant::mathematical::pi*pow3(d_);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::particleShapes::sphere::write
(
    Ostream& os,
    const particleShape& shape
) const
{
    const sphere& s = static_cast<const particleShapes::sphere&>(shape);
    os  << "d" << token::TAB << s.d_ << token::END_STATEMENT << endl;

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const particleShapes::sphere&)"
    );

    return os;
}
// ************************************************************************* //
