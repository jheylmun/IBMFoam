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
    Foam::particleIBM

Description
    abstract class for a particle used with IBM method

SourceFiles
    particleIBM.C

\*---------------------------------------------------------------------------*/

#ifndef particleIBM_H
#define particleIBM_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOstream.H"
#include "vectorList.H"
#include "symmTensorList.H"
#include "labelVector.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "runTimeSelectionTables.H"
#include "particle.H"
#include "particleShape.H"
#include "interpolation.H"


namespace Foam
{

class particleIBM;

Ostream& operator<<(Ostream&, const particleIBM&);


class particleIBM
:
    public particle
{
// Private data

    //- Size in bytes of the fields
    static const std::size_t sizeofFields_;


protected:
// Protected Data

    //- Reference to mesh
    const polyMesh& mesh_;

    //- Dictionary used to constuct particle
    const dictionary& dict_;

    //- Particle name
    label index_;

    //- Is the partice active
    bool active_;

    //- Particle shape containing mesh
    autoPtr<particleShape> shape_;

    //- Linear velocity of particle
    vector v_;

    //- Angular velocity of particle
    vector omega_;

    //- Particle density
    scalar rho_;

    //- Age of particle
    scalar age_;

    //- Total force acting on the particle
    vector integratedForce_;

    //- Pointer to same particle, but across cyclic or processor
    particleIBM* particlePtr_;


    // Protected member fuctions

        //- Interpolate a given field from the global mesh on the the shellMesh
        template<class type>
        void interpolateFromMesh
        (
            const GeometricField<type, fvsPatchField, surfaceMesh>& fieldF,
            Field<type>& field
        ) const;

        //- Interpolate F field to global mesh given U and rho
        void interpolateFToMesh
        (
            const vectorList& U,
            const vectorList& Uold,
            const vectorList& S,
            volVectorField& F
        ) const;


public:

    TypeName("particleIBM");

    // Constructors

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        particleIBM
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const label& index
        );

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        particleIBM
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const label& index,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        particleIBM
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const label& index,
            const vector& position,
            const label celli
        );

        //- Construct from Istream
        particleIBM
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        particleIBM(const particleIBM& p);

        //- Construct as a copy
        particleIBM(const particleIBM& p, const polyMesh& mesh);


    //- Destructor
    virtual ~particleIBM();


    // Access

        //- Particle name
        inline label index() const;

        //- Is any part of the particle on the mesh
        inline particleShape::locationType onMesh() const;

        //- Is the particle active?
        inline bool active() const;

        //- Center of particle
        //  Different from position of particle class so center of particle can
        //  be located off of mesh
        inline const vector& center() const;

        inline vector& center();

        //- Return the particle mass
        inline scalar r(const vector& x2) const;

        //- return effective diameter
        inline scalar d() const;

        //- Return the particle mass
        inline scalar mass() const;

        //- Const access to linear velocity
        inline const vector& v() const;

        //- Non-constant access to linear velocity
        inline vector& v();

        //- Const access to angular velocity
        inline const vector& omega() const;

        //- Non-constant access to angular velocity
        inline vector& omega();

        //- Copy of particle density
        inline const scalar& rho() const;

        //- Non-constant access to density
        inline scalar& rho();

        //- Return total force acting on particle
        inline const vector& F() const;

        //- Return total force acting on particle
        inline vector& F();

        //- Return the age of the particle
        inline const scalar& age() const;

        //- Return the age of the particle
        inline scalar& age();

        //- Return the pointer to another copy of the particle
        inline particleIBM* particlePtr();


    //- Mpi data access

        //- Center processor
        inline label centerProc() const;

        //- List of shared processors const;
        inline const boolList& neiProcs() const;


    // Edit

        //- Calcluate needed forcing on fluid
        void forcing
        (
            const surfaceVectorField& Uf,
            const surfaceVectorField& Ufold,
            const surfaceVectorField& Sf,
            volVectorField& F
        ) const;

        //- Calculate total stress on particle surface
        void integrateSurfaceStress
        (
            const surfaceSymmTensorField& tauf,
            const surfaceScalarField& p
        );

        //- Move particles in space
        void solve(const scalar& dt);

        //- Compute wall hits
        void wallHit(const fvMesh& mesh);

        //- Calculate and return the drag coefficient
        scalar Cd(const scalar& rhoRef, const vector& UInf) const;

        //- Update weights, cell lists, and mesh
        void update();


    //- Read/Write functions

        // Write fields
        static void writeFields(const Cloud<particleIBM>& c);

        // Read fields
        static void readFields(Cloud<particleIBM>& c);


    // Ostream Operator

        // Write shape
        Ostream& writeShape(Ostream& os) const;

        friend Ostream& operator<<
        (
            Ostream&,
            const particleIBM&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End Foam namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "particleIBMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include    "particleIBMTmp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

