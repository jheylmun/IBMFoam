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

#include "particleShape.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(particleShape, 0);
    defineRunTimeSelectionTable(particleShape, dictionary);
    defineRunTimeSelectionTable(particleShape, copy);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::particleShape::setWeights()
{
    // Total number of global mesh faces
    label nFaces = pMesh_.nInternalFaces();

    label nPtsOnMesh = 0;
    forAll(mesh_, pti)
    {
        vector pt = mesh_[pti];
        label celli = pMesh_.findCell(pt);

        if (celli != -1)
        {
            const labelList& cellFaces = pMesh_.cells()[celli];

            //- Local copies of weights and faces
            scalar W = 0.0;
            scalarList w(20, 0.0);
            labelList faces(pMesh_.nFaces(), -1);

            label wi = 0;
            forAll(cellFaces, facei)
            {
                if (cellFaces[facei] >= nFaces) continue;

                faces[wi] = cellFaces[facei];

                scalar diff =
                    mag
                    (
                        pMesh_.faceCentres()[cellFaces[facei]] - pt
                    );

                if (diff < SMALL)
                {
                    w = 0.0;
                    w[wi] = 1.0;
                    W = 1.0;

                    wi++;
                    break;
                }
                w[wi] = (1.0/diff);
                W += w[wi];
                wi++;
            }
            faces.resize(wi);
            w.resize(wi);

            wToLocal_[pti] = w;
            WToLocal_[pti] = W;
            facesToLocal_[pti] = faces;

            nPtsOnMesh++;
        }
    }

    //- Set weighting for interpolation from local mesh
    wFromLocal_ = List<scalarList>(shellCells_.size());
    WFromLocal_ = scalarList(shellCells_.size());

    if (nPtsOnMesh == 0 && pMesh_.findCell(center_) == -1)
    {
        onMesh_ = OFF_MESH;
    }
    else if (nPtsOnMesh == mesh_.size())
    {
        onMesh_ = ON_MESH;
    }
    else
    {
        onMesh_ = PARTIAL;
    }

    forAll(shellCells_, celli)
    {
        label i = shellCells_[celli];
        vector CC = pMesh_.cellCentres()[i];

        scalar W = 0.0;
        scalarList w(8,0.0);

        forAll(w, pti)
        {
            scalar diff = Foam::mag(CC - mesh_[Is_[celli][pti]]);
            if (diff < SMALL)
            {
                w = 0.0;
                w[pti] = 1.0;
                W = 1.0;
                break;
            }
            w[pti] = 1.0/diff;
            W += w[pti];
        }
        wFromLocal_[celli] = w;
        WFromLocal_[celli] = W;
    }
}

void Foam::particleShape::setNeighbours()
{
    Is_ = List<labelList>(shellCells_.size(), labelList(8, 0));
    Os_ = List<labelList>(shellCells_.size(), labelList(4, 0));
    forAll(shellCells_, celli)
    {
        label j = neighbourPoints_[celli].y();
        label k = neighbourPoints_[celli].z();
        label jp1 = ((j+1) % nTheta_);

        if (nk_ != 1)
        {
            if (k == nk_ - 1) continue;

            Is_[celli][0] = index(0, j,   k);
            Is_[celli][1] = index(0, jp1, k);
            Is_[celli][2] = index(0, j,   k + 1);
            Is_[celli][3] = index(0, jp1, k + 1);

            Is_[celli][4] = index(1, j,   k);
            Is_[celli][5] = index(1, jp1, k);
            Is_[celli][6] = index(1, j,   k + 1);
            Is_[celli][7] = index(1, jp1, k + 1);

            Os_[celli][0] = index(2, j,   k);
            Os_[celli][1] = index(2, jp1, k);
            Os_[celli][2] = index(2, j,   k + 1);
            Os_[celli][3] = index(2, jp1, k + 1);
        }
        else
        {
            Is_[celli][0] = index(0, j,   k);
            Is_[celli][1] = index(0, jp1, k);
            Is_[celli][2] = index(0, j,   k);
            Is_[celli][3] = index(0, jp1, k);

            Is_[celli][4] = index(1, j,   k);
            Is_[celli][5] = index(1, jp1, k);
            Is_[celli][6] = index(1, j,   k);
            Is_[celli][7] = index(1, jp1, k);

            Os_[celli][0] = index(2, j,   k);
            Os_[celli][1] = index(2, jp1, k);
            Os_[celli][2] = index(2, j,   k);
            Os_[celli][3] = index(2, jp1, k);
        }
    }
}


void Foam::particleShape::setRotationMatrix()
{
    rotationMatrix_ = tensor::I;

    if (mag(theta_.x()) > SMALL)
    {
        scalar t1 = theta_.x();
        rotationMatrix_ &=
            tensor
            (
                1.0, 0.0, 0.0,
                0.0, cos(t1), -sin(t1),
                0.0, sin(t1), cos(t1)
            );
    }
    if (mag(theta_.y()) > SMALL)
    {
        scalar t2 = theta_.y();
        rotationMatrix_ &=
            tensor
            (
                cos(t2), 0.0, sin(t2),
                0.0, 1.0, 0.0,
                -sin(t2), 0.0, cos(t2)
            );
    }
    if (mag(theta_.z()) > SMALL)
    {
        scalar t3 = theta_.z();
        rotationMatrix_ &=
            tensor
            (
                cos(t3), -sin(t3), 0.0,
                sin(t3), cos(t3), 0.0,
                0.0, 0.0, 1.0
            );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleShape::particleShape
(
    const polyMesh& mesh,
    const dictionary& dict,
    const vector& center,
    const bool buildMesh
)
:
    pMesh_(mesh),
    buildMesh_(buildMesh),
    onMesh_(ON_MESH),
    center_(center),
    momentOfInertia_(HUGE),
    nTheta_(readLabel(dict.lookup("nTheta"))),
    nk_(dict.lookupOrDefault<label>("nk", 1)),
    N_(nRadial_*nTheta_*nk_),
    theta_(dict.lookupOrDefault<vector>("theta", Zero)),
    delta_(readScalar(dict.lookup("delta"))),
    centeredMesh_(N_, Zero),
    mesh_(N_,Zero),
    Sf_(nTheta_*nk_,Zero),
    neighbourPoints_(N_),
    wToLocal_(N_),
    WToLocal_(N_),
    facesToLocal_(N_)
{}


Foam::particleShape::particleShape
(
    const particleShape& shape,
    const vector& center,
    const vector& theta,
    const bool buildMesh
)
:
    pMesh_(shape.pMesh_),
    buildMesh_(buildMesh),
    onMesh_(PARTIAL),
    center_(center),
    momentOfInertia_(shape.momentOfInertia_),
    nTheta_(shape.nTheta_),
    nk_(shape.nk_),
    N_(shape.N_),
    theta_(theta),
    delta_(shape.delta_),
    centeredMesh_(shape.centeredMesh_),
    mesh_(N_, Zero),
    Sf_(shape.Sf_),
    neighbourPoints_(N_),
    wToLocal_(N_),
    WToLocal_(N_),
    facesToLocal_(N_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleShape::~particleShape()
{}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void Foam::particleShape::moveMesh(const vector& center)
{
    center_ = center;
    setRotationMatrix();

    mesh_ = rotate();

    forAll(mesh_, celli)
    {
        mesh_[celli] += center_;
    }
}

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const particleShape& shape
)
{
    return shape.write(os, shape);
}

// ************************************************************************* //
