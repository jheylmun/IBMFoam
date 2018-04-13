/* ---------------------------------------------------------------------------*\
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

#include "gammaDistribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace randomDistributions
{
    defineTypeNameAndDebug(gamma, 0);

    addToRunTimeSelectionTable
    (
        randomDistribution,
        gamma,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::randomDistributions::gamma::gamma
(
    const dictionary& dict
)
:
    randomDistribution(dict),
    alpha_(readScalar(dict.lookup("alpha")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::randomDistributions::gamma::~gamma()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::randomDistributions::gamma::RV(Random& rv)
{
    return invIncGamma(alpha_, rv.scalar01());
}

Foam::scalar Foam::randomDistributions::gamma::moment(const label i) const
{
    return 1.0;
}

Foam::scalar Foam::randomDistributions::gamma::moment(const scalar i) const
{
    return 1.0;
}

// ************************************************************************* //
