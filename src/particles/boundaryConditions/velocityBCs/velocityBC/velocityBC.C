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

#include "velocityBC.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityBC, 0);
    defineRunTimeSelectionTable(velocityBC, dictionary);
    defineRunTimeSelectionTable(velocityBC, copy);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityBC::velocityBC
(
    const dictionary& dict,
    const vectorList& mesh,
    const labelList& shellCells,
    const List<labelList>& Is,
    const List<labelList>& Os,
    const List<scalarList>& ws,
    const scalarList& Ws
)
:
    mesh_(mesh),
    shellCells_(shellCells),
    Is_(Is),
    Os_(Os),
    ws_(ws),
    Ws_(Ws)
{}


Foam::velocityBC::velocityBC(const velocityBC& bc)
:
    mesh_(bc.mesh_),
    shellCells_(bc.shellCells_),
    Is_(bc.Is_),
    Os_(bc.Os_),
    ws_(bc.ws_),
    Ws_(bc.Ws_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityBC::~velocityBC()
{}


// ************************************************************************* //
