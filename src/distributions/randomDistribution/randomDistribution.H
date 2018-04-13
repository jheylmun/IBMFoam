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
    Foam::randomDistribtion

Description
    abstract class random distribtion

SourceFiles
    randomDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef randomDistribution_H
#define randomDistribution_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "Random.H"


namespace Foam
{

class randomDistribution
{
    // Private member functions

        //- Disallow default bitwise copy construct
        randomDistribution(const randomDistribution&);

        //- Disallow default bitwise assignment
        void operator=(const randomDistribution&);


protected:
// Protected Data

    //- Dictionary
    const dictionary& dict_;

    //- Minimum value
    scalar minVal_;

    //- Maximum value
    scalar maxVal_;


public:

    //- Runtime type information
    TypeName("randomDistribution");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        randomDistribution,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    // Selector
    static autoPtr<randomDistribution> New
    (
        const dictionary& dict
    );


    // Constructor
    randomDistribution(const dictionary& dict);


    //- Destructor
    virtual ~randomDistribution();


    // Public functions

        //- Return a random variable based of the given distribution
        virtual scalar RV(Random& rv) = 0;

        //- Return the ith statistical moment
        virtual scalar moment(const label i) const = 0;

        //- Return the ith non-integer moment
        virtual scalar moment(const scalar i) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End Foam namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
