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

const Foam::vectorList& Foam::IBM::particleShape::mesh() const
{
    return baseMesh_;
}

const Foam::vectorList& Foam::IBM::particleShape::Sf() const
{
    return Sf_;
}

Foam::scalar Foam::IBM::particleShape::delta() const
{
    return delta_;
}

Foam::label Foam::IBM::particleShape::N() const
{
    return N_;
}

Foam::label Foam::IBM::particleShape::nRadial() const
{
    return nRadial_;
}

Foam::label Foam::IBM::particleShape::nTheta() const
{
    return nTheta_;
}

Foam::label Foam::IBM::particleShape::nk() const
{
    return nk_;
}

Foam::label Foam::IBM::particleShape::centerIndex() const
{
    return centerIndex_;
}

const Foam::labelList& Foam::IBM::particleShape::shellCells() const
{
    return shellCells_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //