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


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

Foam::particleShape::locationType Foam::particleShape::onMesh() const
{
    return onMesh_;
}

const Foam::vectorList& Foam::particleShape::mesh() const
{
    return mesh_;
}

Foam::tmp<Foam::vectorField> Foam::particleShape::Sf() const
{
    return rotationMatrix_ & Sf_;
}

const Foam::scalar& Foam::particleShape::delta() const
{
    return delta_;
}

Foam::label Foam::particleShape::N() const
{
    return N_;
}

Foam::label Foam::particleShape::nRadial() const
{
    return nRadial_;
}

Foam::label Foam::particleShape::nTheta() const
{
    return nTheta_;
}

Foam::label Foam::particleShape::nk() const
{
    return nk_;
}

const Foam::labelList& Foam::particleShape::shellCells() const
{
    return shellCells_;
}

const Foam::vector& Foam::particleShape::theta() const
{
    return theta_;
}

Foam::vector& Foam::particleShape::theta()
{
    return theta_;
}

const Foam::scalar& Foam::particleShape::I() const
{
    return momentOfInertia_;
}

const Foam::vector& Foam::particleShape::center() const
{
    return center_;
}

Foam::vector& Foam::particleShape::center()
{
    return center_;
}

const Foam::tensor& Foam::particleShape::rotationMatrix() const
{
    return rotationMatrix_;
}

// ************************************************************************* //
