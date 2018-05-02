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
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::velocityBC> Foam::velocityBC::New
(
    const dictionary& dict,
    const vectorList& mesh,
    const labelList& shellCells,
    const List<labelList>& Is,
    const List<labelList>& Os,
    const List<scalarList>& ws,
    const scalarList& Ws
)
{
    const word modelType(dict.lookup("velocityBC"));

    Info<< "Selecting velocity boundary condition " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown velocityBC type "
            << modelType << nl << nl
            << "Valid velocityBCs are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<velocityBC>
        (cstrIter()(dict, mesh, shellCells, Is, Os, ws, Ws));
}


Foam::autoPtr<Foam::velocityBC> Foam::velocityBC::New(const velocityBC& bc)
{
    const word modelType(bc.type());

    copyConstructorTable::iterator cstrIter =
        copyConstructorTablePtr_->find(modelType);

    if (cstrIter == copyConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown velocityBC type "
            << modelType << nl << nl
            << "Valid velocityBCs are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<velocityBC>
        (cstrIter()(bc));
}


// ************************************************************************* //