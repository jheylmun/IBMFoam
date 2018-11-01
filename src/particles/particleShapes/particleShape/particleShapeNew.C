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
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::particleShape> Foam::particleShape::New
(
    const polyMesh& mesh,
    const dictionary& dict,
    const vector& center,
    const bool buildMesh
)
{
    const word modelType(dict.lookup("particleShape"));

    Info<< "Selecting particle shape " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown particle shape type "
            << modelType << nl << nl
            << "Valid particleShapes are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<particleShape>
        (cstrIter()(mesh, dict, center, buildMesh));
}


Foam::autoPtr<Foam::particleShape> Foam::particleShape::New
(
    const particleShape& shape,
    const vector& center,
    const vector& theta,
    const bool buildMesh
)
{
    const word modelType(shape.type());

    copyConstructorTable::iterator cstrIter =
        copyConstructorTablePtr_->find(modelType);

    if (cstrIter == copyConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown particle shape type "
            << modelType << nl << nl
            << "Valid particleShapes are : " << endl
            << copyConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<particleShape>
        (cstrIter()(shape, center, theta, buildMesh));
}


// ************************************************************************* //
