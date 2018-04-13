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
    Foam::particleSizeDistribution

Description
    abstract class for a particleSizeDistribution

SourceFiles
    particleSizeDistribution.C
    particleSizeDistributionNew.C

\*---------------------------------------------------------------------------*/

#ifndef particleSizeDistribution_H
#define particleSizeDistribution_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "vectorList.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "IOstreams.H"
#include "OFstream.H"
#include "Random.H"
#include "constants.H"
#include "boundBox.H"
#include "particleShape.H"
#include "randomDistribution.H"


namespace Foam
{
namespace IBM
{

class particleSizeDistribution
:
    public IOobject
{
    // Private member functions

        //- Disallow default bitwise copy construct
        particleSizeDistribution(const particleSizeDistribution&);

        //- Disallow default bitwise assignment
        void operator=(const particleSizeDistribution&);


protected:
// Protected Data

    //- fvMesh
    const fvMesh& mesh_;

    //- Random number generator
    Random rand_;

    //- Initialization dictionary
    const dictionary& dict_;

    //- List of particle shapes
    PtrList<particleShape> pShapes_;

    //- Types of particles
    wordList pTypes_;

    //- Total number of particles
    label nTotParticles_;

    //- Number of particles in each type
    labelList nParticles_;

    //- Bounding box
    boundBox bb_;

    //- List of particle densities
    scalarField rhos_;


public:

    //- Runtime type information
    TypeName("particleSizeDistribution");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        particleSizeDistribution,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (mesh, dict)
    );

    // Selectors

        static autoPtr<particleSizeDistribution> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        // Construct particle from dictionary entry
        particleSizeDistribution
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~particleSizeDistribution();


    // IO functions

        bool write();

        virtual bool writeParticles(Ostream& os) = 0;
};

}   // End IBM namespace
}   // End Foam namespace


#endif
