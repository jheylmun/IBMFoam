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
namespace IBM
{
    defineTypeNameAndDebug(particleShape, 0);
    defineRunTimeSelectionTable(particleShape, dictionary);
}
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::IBM::particleShape::setWeights()
{
    wFromLocal_.clear();

    // Total number of global mesh faces
    label nFaces = mesh_.nInternalFaces();

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

                if (diff < small)
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
        }
    }

    //- Set weighting for interpolation from local mesh
    wFromLocal_ = List<scalarList>(shellCells_.size());
    WFromLocal_ = scalarList(shellCells_.size());

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

    centerIndex_ = mesh_.findCell(position());
}

void Foam::IBM::particleShape::setNeighbours()
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IBM::particleShape::particleShape
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    nTheta_(readLabel(dict.lookup("nTheta"))),
    nk_(dict.lookupOrDefault<label>("nk",1)),
    N_(nRadial_*nTheta_*nk_),
    delta_(readScalar(dict.lookup("delta"))),
    centeredMesh_(N_, Zero),
    baseMesh_(N_,Zero),
    Sf_(nTheta_*nk_,Zero),
    neighbourPoints_(N_),
    wToLocal_(N_),
    WToLocal_(N_),
    facesToLocal_(N_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IBM::particleShape::~particleShape()
{}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void Foam::IBM::particleShape::moveMesh()
{
    forAll(baseMesh_, celli)
    {
        baseMesh_[celli] = centeredMesh_[celli] + position();
    }
}