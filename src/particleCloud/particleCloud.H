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
    Foam::particleCloud

Description
    Templated list of particles

SourceFiles
    particleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef particleCloud_H
#define particleCloud_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "vectorList.H"
#include "labelVector.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "runTimeSelectionTables.H"
#include "particleIBM.H"
#include "Cloud.H"
#include "wallDist.H"


namespace Foam
{

template<class cloudType>
class IOposition;

class particleCloud;

Ostream& operator<<
(
    Ostream&,
    const particleCloud&
);


class particleCloud
:
    public Cloud<particleIBM>
{
    // Private Data

        // Constant access to fvMesh
        const fvMesh& mesh_;

        //- dictionary used to initialize points
        dictionary dict_;

        //- Number of particles since simuation start
        label nParticles_;

        //- List of original particles
        UPtrList<particleIBM> origParticles_;

        //- Reference values
        vector UInf_;
        vector flowNormal_;
        scalar rhoRef_;
        scalar pRef_;

        //- Gravitational acceleration
        vector g_;

        //- Coefficient of restitution
        scalar e_;

        //- Copy of old velocity field
        volVectorField Uold_;

        //- IBM source terms (advection, pressure, and viscous stress)
        volVectorField S_;

        //- Is the solution coupled to fluid
        Switch coupled_;

        //- Switch if particles are moving
        Switch moving_;

        //- Switch if particles are rotating
        Switch rotating_;

        //- List of patch numbers used by particles
        //  walls, cyclic, processor, and cyclicProcessor
        labelList patches_;

        //- Are there walls on the mesh
        labelList walls_;

        //- Maximum copy number
        labelList maxCopies_;

        //- Indicator function
        volScalarField I_;

        // Collision variables

    // Private functions

        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volSymmTensorField> devRhoReff() const;

        //- Dynamic viscosity field
        tmp<volScalarField> mu() const;

        //- Return rho if specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Compute forces due to particle collisions
        void computeCollisions(const scalar& dt);

        //- Compute forces due to wall collisions
        void computeWallHits(const scalar& dt);

        //- Does a particle already exist
        bool found
        (
            const label index,
            const label copy,
            const vector& center
        ) const;

        //- Return const particle pointer if it exists
        const particleIBM* find(const label index, const label copy) const;

        //- Return non-const particle pointer if it exists
        particleIBM* find(const label index, const label copy);

        //- Check if particles cross cyclic boundaries
        void computeCyclicHits();

        //- Set sorces of particles that have points on mesh,
        void combineForces();

        //- Set positions of particles that have points on mesh,
        //  but not the center
        void setPositions();

        //- Delete particles no longer on mesh
        void deleteInActiveParticles();

        //- Update the indicator function field
        void updateI();


    // Combine fields ops

        //- Add over procs
        template<class T>
        class plusEqOp
        {

        public:

            void operator()(T& x, const T& y) const
            {
                x += y;
            }
        };

        //- Add over procs
        template<class T>
        class ListPlusEqOp
        {

        public:

            void operator()(T& x, const T& y) const
            {
                x.setSize(x.size());

                forAll(y, i)
                {
                    x[i] += y[i];
                }
            }
        };

        //- Max value across procs
        template<class T>
        class combineListEqOp
        {

        public:

            void operator()(T& x, const T& y) const
            {
                x.append(y);
            }
        };

        //- Max value across procs
        template<class T>
        class ListMaxEqOp
        {

        public:

            void operator()(T& x, const T& y) const
            {
                forAll(x, i)
                {
                    x[i] = max(x[i], y[i]);
                }
            }
        };
        //- Max value across procs
        template<class T>
        class ListListMaxVectorEqOp
        {

        public:

            void operator()(T& x, const T& y) const
            {
                forAll(x, i)
                {
                    forAll(x[i], j)
                    {
                        if (mag(y[i][j]) > mag(x[i][j]))
                        {
                            x[i][j] = y[i][j];
                        }
                    }
                }
            }
        };


public:

    TypeName("particleCloud");
    // Constructors

        // Construct particleCloud from dictionary entry
        particleCloud
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionedVector& g
        );


    //- Destructor
    virtual ~particleCloud();


    // Access

        //- Number of particle
        label nParticles() const
        {
            return refCast<const Cloud<particleIBM>>(*this).size();
        }

        //- Const access to fvMesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return dictionary
        const dictionary& ibmDict() const
        {
            return dict_;
        }

        //- Freestream velocity
        vector UInf() const
        {
            return UInf_;
        }

        //- Freestream velocity
        scalar rhoRef() const
        {
            return rhoRef_;
        }

        //- Freestream velocity
        scalar pRef() const
        {
            return pRef_;
        }

        //- Return particle courant number (v*deltaT/delta)
        scalar maxCoNum() const;

    // Edit

        //- Calcluate needed forcing on fluid
        tmp<volVectorField> forcing
        (
            const volVectorField& S
        );

        //- Integrate stresses on particle surfaces
        void integrateSurfaceStresses();

        //- Calculate drag coefficients
        void Cd(IOField<scalar>& Cds) const;

        //- Advance in time, update velocities based on force from fluid, and
        //  move particles in space.
        void solve();


    // Read


    //- Write functions

        //- Write the field data for the cloud of particles Dummy at
        //  this level.
        virtual void writeFields() const;

        //- Write using given format, version and compression.
        //  Only writes the cloud file if the Cloud isn't empty
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool valid
        ) const;

        void readParticles();


    // Ostream operator

        friend Ostream& operator<<
        (
            Ostream&,
            const particleCloud&
        );


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End Foam namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
