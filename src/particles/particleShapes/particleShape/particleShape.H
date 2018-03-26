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

#ifndef particleShape_H
#define particleShape_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "vectorList.H"
#include "labelVector.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "polyMesh.H"
#include "constants.H"
#include "IOstream.H"


namespace Foam
{

class particleIBM;
class particleShape;

Ostream& operator<<(Ostream&, const particleShape&);

class particleShape
{
public:

    //- Declaration of on mesh indicator
    enum locationType
    {
        OFF_MESH,
        PARTIAL,
        ON_MESH
    };


private:
    // Private member functions

        //- Disallow default bitwise copy construct
        particleShape(const particleShape&);

        //- Disallow default bitwise assignment
        void operator=(const particleShape&);


protected:
// Protected Data

    friend particleIBM;

    //- fvMesh
    const polyMesh& mesh_;

    //- dictionary
    const dictionary& dict_;

    //- Is the particle on the mesh
    locationType onMesh_;

    //- center point
    vector center_;

    //- number of points in the radial direction (3)
    label nRadial_ = 3;

    //- Number of points in the theta direction
    label nTheta_;

    //- Number of points in third direction (phi or z)
    label nk_;

    //- Total number of grid points
    label N_;

    //- Delta that is used as the spacing from the actual partical diameter,
    //  ie x_o = pos + delta
    scalar delta_;

    //- Mesh centered at (0,0,0)
    vectorList centeredMesh_;

    //- Base mesh since particles do not deform
    vectorList baseMesh_;

    //- Surface normal vector time surface area
    vectorList Sf_;

    //- Cell containing center of shape
    label centerIndex_;

    //- List of cells within the particles shells
    mutable labelList shellCells_;

    //- Contains the lowest integer that surrounds the cell center
    List<labelVector> neighbourPoints_;

    //- List containing the indicies of the neighbour points
    //  (inner and surface)
    List<labelList> Is_;

    //- List containing the indicies of the neighbour points
    //  (outer)
    List<labelList> Os_;

    //- Weighting parameters for interpolation from global mesh
    mutable List<scalarList> wToLocal_;

    //- Total weighting parameters for interpolation from global mesh
    mutable scalarList WToLocal_;

     //- Face numbers for interpolation from global mesh
    mutable List<labelList> facesToLocal_;

    //- Weighting parameters for interpolation to global mesh
    mutable List<scalarList> wFromLocal_;

    //- Total weighting parameters for interpolation to global mesh
    mutable scalarList WFromLocal_;


    // Mpi data

        //- processor id at center
        label centerProc_;

        //- Processors used in interpolation of surface mesh
        boolList neiProcs_;


    // Protected member fuctions

        //- Return total index based on theta, phi index
        label index(const int i, const int j, const int k) const
        {
            return
            (
                i
              + j*nRadial_
              + k*nRadial_*nTheta_
            );
        }
        label index2(const int i, const int j) const
        {
            return
            (
                i
              + j*nTheta_
            );
        }

        //- Set weighting values for interpolation
        void setWeights();

        //- Set Is and Os
        void setNeighbours();

        //- Set surface areas and surface normals
        virtual void calcSf() = 0;

        //- Set owner and neighbour processors
        void setProcs();


public:

    //- Runtime type information
    TypeName("particleShape");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        particleShape,
        dictionary,
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const vector& center
        ),
        (mesh, dict, center)
    );

    // Selectors

        static autoPtr<particleShape> New
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const vector& center
        );


    // Constructors

        // Construct particle from dictionary entry
        particleShape
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const vector& center
        );


    //- Destructor
    virtual ~particleShape();


    // Access

        //- Is the shape on the mesh
        inline locationType onMesh() const;

        //- Is the particle on a single processor
        inline bool singleProc() const;

        //- Return structured mesh
        inline const vectorList& mesh() const;

        //- Return surface normals
        inline const vectorList& Sf() const;

        //- Return delta
        inline const scalar& delta() const;

        //- Return total number of points
        inline label N() const;

        //- Return number of points in radial direction
        inline label nRadial() const;

        //- Return number of points in theta direction
        inline label nTheta() const;

        //- Return number of points in k or phi direction
        inline label nk() const;

        //- Return index containing cell at center of shape
        inline label centerIndex() const;

        //- Return cell in the shell of the surface mesh
        inline const labelList& shellCells() const;

        //- Particle diameter, or equivelant diameter
        virtual scalar d() const = 0;

        //- Frontal Area
        virtual scalar A() const = 0;

        //- Particle volume
        virtual scalar V() const = 0;

        //- Return center of mass
        inline const vector& center() const;

        //- Return center of mass
        inline vector& center();

        //- Return the radius give a position of another particle
        virtual scalar r(const vector& x2) const = 0;


    // Edit functions

        //- Update onMesh variable
        void setOnMesh();

        //- Update cells inside particle and shells
        virtual void updateCellLists() = 0;

        //- Discretize particle on the surface and return the list of points.
        //  The created mesh will store values at the verticies rather than the
        //  cell centers.
        virtual void discretize() = 0;

        //- Move mesh from (0,0,0) to CoM
        virtual void moveMesh(const vector& center);


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const particleShape&
        );

        // Write
        virtual Ostream& write
        (
            Ostream&,
            const particleShape&
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End Foam namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include    "particleShapeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
