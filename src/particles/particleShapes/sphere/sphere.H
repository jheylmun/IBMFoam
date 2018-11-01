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

#ifndef sphere_H
#define sphere_H

#include "particleShape.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace particleShapes
{

class sphere
:
    public particleShape
{
    // Private member functions

        //- Particle diameter
        scalar d_;

        //- Point1
        vector center_;

        //- Calculate surface normal vector time area
        virtual void calcSf();


public:

    //- Runtime type information
    TypeName("sphere");

    // Constructor

        // Construct particle from dictionary entry
        sphere
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const vector& center,
            const bool buildMesh = true
        );

        sphere
        (
            const particleShape& shape,
            const vector& center,
            const vector& theta,
            const bool buildMesh = true
        );

    //- Destructor
    virtual ~sphere();


    // Member functions

        //- Update cells inside particle and shells
        virtual void updateCellLists();

        //- Discretize particle on the surface and return the list of points.
        //  The created mesh will store values at the verticies rather than the
        //  cell centers.
        virtual void discretize();

        //- Particle diameter
        virtual scalar d() const;

        //- Tensor of diameters
        virtual vector D() const;

        //- Surface area
        virtual scalar A(const vector& pt) const;

        //- Particle volume
        virtual scalar V() const;

        //- Return the radius give a position of another particle
        virtual scalar r(const vector& x2) const
        {
            return d_/2.0;
        }


    // Ostream Operator

        virtual Ostream& write
        (
            Ostream&,
            const particleShape&
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End particleShapes namespace
}   // End Foam namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
