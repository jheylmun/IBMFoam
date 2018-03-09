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
    Foam::particleList

Description
    Templated list of particles

SourceFiles
    particleList.C

\*---------------------------------------------------------------------------*/

#ifndef particleList_H
#define particleList_H

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
#include "particle.H"
#include "regIOobject.H"
#include "cloud.H"


namespace Foam
{
namespace IBM
{

template<class pType>
class particleList
:
    public regIOobject,
    public PtrList<pType>
{
    // Private Data

        // Constant access to fvMesh
        const fvMesh& mesh_;

        //- dictionary used to initialize points
        IOdictionary dict_;

        //- Reference values
        vector UInf_;
        vector flowNormal_;
        scalar rhoRef_;
        scalar pRef_;

        //- Coefficient of restitution
        scalar e_;

        //- Copy of old velocity field
        volVectorField Uold_;

        //- IBM source terms (advection, pressure, and viscous stress)
        volVectorField S_;

        // Switche if particles are moving
        Switch moving_;


    // Priavte functions

        //- Create IOobject
        IOobject fieldIOobject
        (
            const word& fileName,
            const IOobject::readOption& r
        ) const;

        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volSymmTensorField> devRhoReff() const;

        //- Dynamic viscosity field
        tmp<volScalarField> mu() const;

        //- Return rho if specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Calculate collisions
        void computeCollisions();


public:

    // Constructors

        // Construct particleList from dictionary entry
        particleList
        (
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~particleList();


    // Access

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

        //- Normalized freestream velocity
        vector normUInf() const
        {
            return UInf_/mag(UInf_);
        }

        //- Freestream velocity
        scalar rhoRef() const
        {
            return rhoRef_;
        }

        //- Freestream velocity
        scalar pRef() const
        {
            return pRef;
        }

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
        void operator++();

    //- regIOobject functions

        virtual bool write() const;
        virtual bool writeData(Ostream& os) const;


};

}   // End IBM namespace

}   // End Foam namespace

#if NoRepository
#include    "particleList.C"
#endif

#endif
