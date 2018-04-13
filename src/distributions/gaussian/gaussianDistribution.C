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

#include "gaussianDistribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace randomDistributions
{
    defineTypeNameAndDebug(gaussian, 0);

    addToRunTimeSelectionTable
    (
        randomDistribution,
        gaussian,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::randomDistributions::gaussian::gaussian
(
    const dictionary& dict
)
:
    randomDistribution(dict),
    mean_(readScalar(dict.lookup("mean"))),
    variance_(readScalar(dict.lookup("variance")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::randomDistributions::gaussian::~gaussian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::randomDistributions::gaussian::RV(Random& rv)
{
    return rv.GaussNormal()*variance_ + mean_;
}

Foam::scalar Foam::randomDistributions::gaussian::moment(const label i) const
{
    if (i == 0)
    {
        return 1.0;
    }
    else if (i == 1)
    {
        return mean_;
    }
    else if (i == 2)
    {
        return sqr(mean_) + sqr(variance_);
    }
    else if (i == 3)
    {
        return mean_*(sqr(mean_) + 2.0*sqr(variance_));
    }
    else if (i == 4)
    {
        return pow4(mean_) + 6.0*sqr(mean_)*sqr(variance_) + pow4(variance_);
    }

    NotImplemented;
    return 1.0;
}

Foam::scalar Foam::randomDistributions::gaussian::moment(const scalar i) const
{
    NotImplemented;
    return 1.0;
}

// ************************************************************************* //
