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

#ifndef particle_H
#define particle_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "vectorList.H"
#include "symmTensorList.H"
#include "labelVector.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "runTimeSelectionTables.H"
#include "polyMesh.H"
#include "particleShape.H"
#include "interpolation.H"


namespace Foam
{
namespace IBM
{

template<class pType>
class particleList;


class particle
{
protected:
// Protected Data

    //- Constant access to list of particles
    const particleList<particle>& pList_;

    //- fvMesh
    const fvMesh& mesh_;

    //- Particle shape containing mesh
    autoPtr<particleShape> shape_;

    //- Const ref to velocity
    const volVectorField& U_;

    //- Const ref to old velocity
    const volVectorField& U0_;

    //- Const ref to density
//     const volScalarField& rhof_;

    //- Const ref to pressure
    const volScalarField& p_;

    //- Const ref to IBM source term
    const volVectorField& S_;

//     //- Velocity interpolator
//     autoPtr<interpolation<vector>> UInterp_;
//
//     //- Velocity_old interpolator
//     autoPtr<interpolation<vector>> U0Interp_;
//
//     //- Density interpolator
//     autoPtr<interpolation<scalar>> rhoInterp_;
//
//     //- Pressure interpolator
//     autoPtr<interpolation<scalar>> pInterp_;
//
//     //- Source interpolator
//     autoPtr<interpolation<vector>> SInterp_;

    //- Linear velocity of particle
    vector v_;

    //- Angular velocity of particle
    vector omega_;

    //- Particle density
    scalar rho_;

    //- Total force acting on the particle
    vector integratedForce_;

    // Mpi data

        //- processor id at center
        label origProc_;

        //- Processors used in interpolation of surface mesh
        labelList neiProcs_;


    // Protected member fuctions

        //- Interpolate a given field from the global mesh on the the shellMesh
        template<class type>
        void interpolateFromMesh
        (
            const GeometricField<type, fvsPatchField, surfaceMesh>& fieldF,
            List<type>& field
        ) const;

        //- Interpolate F field to global mesh given U and rho
        void interpolateFToMesh
        (
            const vectorList& U,
            const vectorList& Uold,
            const vectorList& S,
            volVectorField& F
        );


        //- Set owner and neighbour processors
        void setProcs();


public:

    // Constructors

        // Construct particle from dictionary entry
        particle
        (
            const particleList<particle>& pList,
            const dictionary& dict,
            const volVectorField& Uold,
            const volVectorField& S
        );


    //- Destructor
    virtual ~particle();


    // combine fields
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
        template<class T>
        class MaxEqOp
        {

        public:

            void operator()(T& x, const T& y) const
            {
                x = max(x,y);
            }
        };
    // Access

        //- Const access to particle list
        const particleList<particle>& pList() const
        {
            return pList_;
        }

        //- Const access to linear velocity
        vector v() const
        {
            return v_;
        }

        //- Non-constant access to linear velocity
        vector& v()
        {
            return v_;
        }

        //- Const access to angular velocity
        vector omega() const
        {
            return omega_;
        }

        //- Non-constant access to angular velocity
        vector& omega()
        {
            return omega_;
        }

        //- Copy of particle density
        scalar rho() const
        {
            return rho_;
        }

        //- Non-constant access to density
        scalar& rho()
        {
            return rho_;
        }

        const vector& CoM() const
        {
            return shape_->CoM();
        }

        vector& CoM()
        {
            return shape_->CoM();
        }

    // Edit

        //- Calcluate needed forcing on fluid
        void forcing
        (
            const surfaceVectorField& Uf,
            const surfaceVectorField& Ufold,
            const surfaceVectorField& Sf,
            volVectorField& F
        );

        void integrateSurfaceStress
        (
            const surfaceSymmTensorField& tauf,
            const surfaceScalarField& p
        );

        scalar Cd() const;



        //- Advance in time, update velocity based on force from fluid, and
        //  move particle in space.
        virtual void operator++()
        {}

        //- Update weights, cell lists, and mesh
        void update();
};

}   // End IBM namespace

}   // End Foam namespace

#if NoRepository
#include    "particleTmp.H"
#endif

#endif
