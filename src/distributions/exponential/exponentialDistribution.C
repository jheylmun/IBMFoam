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

#include "exponentialDistribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace randomDistributions
{
    defineTypeNameAndDebug(exponential, 0);

    addToRunTimeSelectionTable
    (
        randomDistribution,
        exponential,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::randomDistributions::exponential::exponential
(
    const dictionary& dict
)
:
    randomDistribution(dict),
    mean_(readScalar(dict.lookup("mean")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::randomDistributions::exponential::~exponential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::randomDistributions::exponential::RV(Random& rv)
{
    return -log(rv.scalar01())/mean_;
}

Foam::scalar Foam::randomDistributions::exponential::moment(const label i) const
{
    return pow(-1*mean_, i);
}

Foam::scalar Foam::randomDistributions::exponential::moment(const scalar i) const
{
    NotImplemented;
    return 1.0;
}

// ************************************************************************* //
