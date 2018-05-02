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

#ifndef slipBC_H
#define slipBC_H

#include "velocityBC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
namespace velocityBCs
{

class slipBC
:
    public velocityBC
{
public:

    //- Runtime type information
    TypeName("slip");


    // Constructors

        // Construct shape from dictionary entry
        slipBC
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
        slipBC(const velocityBC& bc);


    //- Destructor
    virtual ~slipBC();


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
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End namespace velocityBCs
}   // End Foam namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
