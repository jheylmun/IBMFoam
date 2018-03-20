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

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class T>
void Foam::particleIBM::interpolateFromMesh
(
    const GeometricField<T, fvsPatchField, surfaceMesh>& fieldF,
    List<T>& field
) const
{
    if (centerProc_ == -1 && neiProcs_.size() == 0)
    {
        return;
    }

    const List<scalarList>& ws = shape_->wToLocal_;
    const scalarList& W = shape_->WToLocal_;
    const List<labelList>& facesList = shape_->facesToLocal_;

    forAll(shape_->baseMesh_, pti)
    {
        const labelList& faces = facesList[pti];
        forAll(faces, facei)
        {
            field[pti] += ws[pti][facei]/W[pti]*fieldF[faces[facei]];
        }
    }
}

// ************************************************************************* //
