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
    Foam::exponential

Description
    Random exponential distribution

SourceFiles
    exponentialDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialDistribution_H
#define exponentialDistribution_H

#include "randomDistribution.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace randomDistributions
{

class exponential
:
    public randomDistribution
{

    //- Mean value
    scalar mean_;


public:

    //- Runtime type information
    TypeName("exponential");

    // Constructor
    exponential(const label seed, const dictionary& dict);

    //- Destructor
    virtual ~exponential();


    // Public functions

        //- Return a exponential random variable
        virtual scalar RV();

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
