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
    Foam::movingParticleIBM

Description
    abstract class for a particle used with IBM method

SourceFiles
    movingParticleIBM.C

\*---------------------------------------------------------------------------*/

#ifndef movingParticleIBM_H
#define movingParticleIBM_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "particleIBM.H"


namespace Foam
{

template<class pType> class movingParticleIBM;

Ostream& operator<<(Ostream&, const movingParticleIBM<pType>&);


class movingParticleIBM
:
    public pType
{
// Private data

    //- Size in bytes of the fields
    static const std::size_t sizeofFields_;


protected:
// Protected Data

    //- Linear velocity of particle
    vector v_;


    // Protected member fuctions

        //- Return the surface velocity given a point
        virtual vector v(const vector& pt) const;


public:

    TypeName("moving");

    // Constructors

        //- Construct from a dictionary, used as initial construction of
        //  particle, i.e. copy 0
        movingParticleIBM
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const label index
        );

        //- Construct from Istream
        movingParticleIBM
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields
        );

        //- Construct as a copy with new center (used with cyclic boundaries)
        movingParticleIBM
        (
            const movingParticleIBM<pType>& p,
            const label copy,
            const vector& center,
            const vector& theta,
            const vector& v,
            const vector& omega,
            const scalar& age
        );

        //- Construct as a copy
        movingParticleIBM(const movingParticleIBM<pType>& p);


    //- Destructor
    virtual ~movingParticleIBM();


    // Access

        //- Const access to linear velocity
        virtual const vector& v() const
        {
            return v_;
        }

        //- Non-constant access to linear velocity
        virtual vector& v()
        {
            return v_;
        }


    //- Edit

        //- Calculate total stress on particle surface
        void integrateSurfaceStress
        (
            const surfaceSymmTensorField& tauf,
            const surfaceScalarField& p
        );

        //- Move particles in space
        virtual void solve(const scalar& dt);


        //- Update weights, cell lists, and mesh
        virtual void update();


    //- Read/Write functions

        // Write fields
        static void writeFields
        (
            const typename Cloud<movingParticleIBM<pType>>& c
        );

        // Read fields
        static void readFields
        (
            typename Cloud<movingParticleIBM<pType>>& c
        );


    // Ostream Operator

        // Write shape
        Ostream& writeShape(Ostream& os) const;

        friend Ostream& operator<<
        (
            Ostream&,
            const movingParticleIBM<pType>&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End Foam namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include    "movingParticleIBM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

