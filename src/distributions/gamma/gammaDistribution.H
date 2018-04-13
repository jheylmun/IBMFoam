/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Class
    Foam::gamma

Description
    Random gamma distribution

SourceFiles
    gammaDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef gammaDistribution_H
#define gammaDistribution_H

#include "randomDistribution.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace randomDistributions
{

class gamma
:
    public randomDistribution
{

    //- Scalae parameter
    scalar alpha_;


public:

    //- Runtime type information
    TypeName("gamma");

    // Constructor
    gamma(const dictionary& dict);

    //- Destructor
    virtual ~gamma();


    // Public functions

        //- Return a gamma random variable
        virtual scalar RV(Random& rv);

        //- Return the ith statistical moment
        virtual scalar moment(const label i) const;

        //- Return the ith non-integer moment
        virtual scalar moment(const scalar i) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End randomDistributions namespace
}   // End Foam namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
