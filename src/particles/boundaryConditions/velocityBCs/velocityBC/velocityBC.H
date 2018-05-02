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
    Foam::particle

Description
    abstract class for a spherical particle

SourceFiles
    particle.C

\*---------------------------------------------------------------------------*/

#ifndef velocityBC_H
#define velocityBC_H

#include "particleIBM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class velocityBC
{
protected:
// Protected Data

    //- Centered mesh of particle
    const vectorList& mesh_;

    //- List of cells inside inside shell
    const labelList shellCells_;

    //- List containing the indicies of the neighbour points
    //  (inner and surface)
    const List<labelList>& Is_;

    //- List containing the indicies of the neighbour points
    //  (outer)
    const List<labelList>& Os_;

    //- Weighting parameters for interpolation to global mesh
    const List<scalarList>& ws_;

    //- Total weighting parameters for interpolation to global mesh
    const scalarList& Ws_;



public:

    //- Runtime type information
    TypeName("velocityBC");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        velocityBC,
        dictionary,
        (
            const dictionary& dict,
            const vectorList& mesh,
            const labelList& shellCells,
            const List<labelList>& Is,
            const List<labelList>& Os,
            const List<scalarList>& ws,
            const scalarList& Ws
        ),
        (dict, mesh, shellCells, Is, Os, ws, Ws)
    );

    declareRunTimeSelectionTable
    (
        autoPtr,
        velocityBC,
        copy,
        (
            const velocityBC& bc
        ),
        (bc)
    );

    // Selectors

        static autoPtr<velocityBC> New
        (
            const dictionary& dict,
            const vectorList& mesh,
            const labelList& shellCells,
            const List<labelList>& Is,
            const List<labelList>& Os,
            const List<scalarList>& ws,
            const scalarList& Ws
        );

        static autoPtr<velocityBC> New(const velocityBC& bc);


    // Constructors

        // Construct shape from dictionary entry
        velocityBC
        (
            const dictionary& dict,
            const vectorList& mesh,
            const labelList& shellCells,
            const List<labelList>& Is,
            const List<labelList>& Os,
            const List<scalarList>& ws,
            const scalarList& Ws
        );

        // Construct as a copy
        velocityBC(const velocityBC& bc);


    //- Destructor
    virtual ~velocityBC();


    // Public member functions

        //- Calculate forcing on mesh
        virtual void forcingOnFluid
        (
            const scalar& dT,
            const vectorList& U,
            const vectorList& Uold,
            const vectorList& S,
            const particleIBM& p,
            volVectorField& F
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End Foam namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
