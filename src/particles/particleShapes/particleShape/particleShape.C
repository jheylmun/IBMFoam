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
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::particleShape::setWeights()
{
    wFromLocal_.clear();

    // Total number of global mesh faces
    label nFaces = mesh_.nInternalFaces();

    label nPtsOnMesh = 0;
    forAll(baseMesh_, pti)
    {
        vector pt = baseMesh_[pti];
        label celli = mesh_.findCell(pt);

        if (celli != -1)
        {
            const labelList& cellFaces = mesh_.cells()[celli];

            //- Local copies of weights and faces
            scalar W = 0.0;
            scalarList w(20, 0.0);
            labelList faces(mesh_.nFaces(), -1);

            label wi = 0;
            forAll(cellFaces, facei)
            {
                if (cellFaces[facei] > nFaces) continue;

                faces[wi] = cellFaces[facei];

                scalar diff =
                    mag
                    (
                        mesh_.faceCentres()[cellFaces[facei]] - pt
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

    centerIndex_ = mesh_.findCell(center_);
    if (nPtsOnMesh == 0 && centerIndex_ == -1)
    {
        onMesh_ = OFF_MESH;
    }
    else if (nPtsOnMesh == baseMesh_.size() && centerIndex_ != -1)
    {
        onMesh_ = ON_MESH;
    }
    else if (centerIndex_ != -1)
    {
        onMesh_ = PARTIAL_WITH_CENTER;
    }
    else
    {
        onMesh_ = PARTIAL_WITHOUT_CENTER;
    }

    forAll(shellCells_, celli)
    {
        label i = shellCells_[celli];
        vector CC = mesh_.cellCentres()[i];

        scalar W = 0.0;
        scalarList w(8,0.0);

        forAll(w, pti)
        {
            scalar diff = Foam::mag(CC - baseMesh_[Is_[celli][pti]]);
            if (diff < 1e-8)
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


void Foam::particleShape::setProcs()
{
    if (centerIndex_ != -1)
    {
        centerProc_ = Pstream::myProcNo();
    }

    neiProcs_ = false;
    forAll(shellCells_, celli)
    {
        if
        (
            shellCells_[celli] != -1
         && Pstream::myProcNo() != centerProc_
        )
        {
            neiProcs_[Pstream::myProcNo()] = true;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleShape::particleShape
(
    const polyMesh& mesh,
    const dictionary& dict,
    const vector& center
)
:
    mesh_(mesh),
    dict_(dict),
    onMesh_(ON_MESH),
    center_(center),
    momentOfInertia_(HUGE),
    nTheta_(readLabel(dict.lookup("nTheta"))),
    nk_(dict.lookupOrDefault<label>("nk", 1)),
    N_(nRadial_*nTheta_*nk_),
    theta_(dict.lookupOrDefault<vector>("nk", Zero)),
    delta_(readScalar(dict.lookup("delta"))),
    centeredMesh_(N_, Zero),
    baseMesh_(N_,Zero),
    Sf_(nTheta_*nk_,Zero),
    neighbourPoints_(N_),
    wToLocal_(N_),
    WToLocal_(N_),
    facesToLocal_(N_),
    neiProcs_(Pstream::nProcs(), false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleShape::~particleShape()
{}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void Foam::particleShape::moveMesh(const vector& center)
{
    center_ = center;

    forAll(baseMesh_, celli)
    {
        baseMesh_[celli] = centeredMesh_[celli] + center_;
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
