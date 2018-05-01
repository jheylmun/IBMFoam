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
    Foam::nonRotatingParticleIBM

Description
    abstract class for a particle used with IBM method

SourceFiles
    nonRotatingParticleIBM.C

\*---------------------------------------------------------------------------*/

#ifndef nonRotatingParticleIBM_H
#define nonRotatingParticleIBM_H

#include "particleIBM.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class pType> class nonRotatingParticleIBM;

Ostream& operator<<(Ostream&, const nonRotatingParticleIBM<pType>&);


class nonRotatingParticleIBM
:
    public pType
{
// Private data

    //- Size in bytes of the fields
    static const std::size_t sizeofFields_;


protected:
// Protected Data

    //- Angular velocity of particle
    vector omega_ = Zero;


    // Protected member fuctions

        //- Return the surface velocity given a point from the centered mesh
        virtual vector v(const vector& pt) const;


public:

    TypeName("nonRotatingParticleIBM");

    // Constructors

        //- Construct from a dictionary, used as initial construction of
        //  particle, i.e. copy 0
        nonRotatingParticleIBM
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const label index
        );

        //- Construct from Istream
        nonRotatingParticleIBM
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields
        );

        //- Construct as a copy with new center (used with cyclic boundaries)
        nonRotatingParticleIBM
        (
            const nonRotatingParticleIBM<pType>& p,
            const label copy,
            const vector& center,
            const vector& theta,
            const vector& v,
            const vector& omega,
            const scalar& age
        );

        //- Construct as a copy
        nonRotatingParticleIBM(const nonRotatingParticleIBM<pType>& p);


    //- Destructor
    virtual ~nonRotatingParticleIBM();


    // Access

        //- Const access to angular velocity
        virtual const vector& omega() const
        {
            return omega_;
        }

        //- Non-constant access to angular velocity
        virtual vector& omega()
        {
            return omega_;
        }


    // Public functions

        //- Move particles in space
        virtual void solve(const scalar& dt);

        //- Update weights, cell lists, and mesh
        void update();


    //- Read/Write functions

        // Write fields
        static void writeFields
        (
            const typename Cloud<nonRotatingParticleIBM<pType>>& c
        );

        // Read fields
        static void readFields
        (
            typename Cloud<nonRotatingParticleIBM<pType>>& c
        );


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const nonRotatingParticleIBM<pType>&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End Foam namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include    "nonRotatingParticleIBM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

