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
    monodisperse.C

\*---------------------------------------------------------------------------*/

#ifndef monodisperse_H
#define monodisperse_H

#include "particleSizeDistribution.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace IBM
{

class monodisperse
:
    public particleSizeDistribution
{
    // Private member functions

        //- Particle diameter
        scalar d_;

        //- Bounding box
        vector bb1_;
        vector bb2_;

        //- points for construction of particle
        vector p1_;
        vector p2_;

        //- List of particle positions
        vectorList particlePos_;


public:

    //- Runtime type information
    TypeName("monodisperse");

    // Constructor

        // Construct particle from dictionary entry
        monodisperse
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~monodisperse();


    // IO functions
    virtual bool writeParticles(Ostream& os);

};

}   // End IBM namespace
}   // End Foam namespace

#endif
