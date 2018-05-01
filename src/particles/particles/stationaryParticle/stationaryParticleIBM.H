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
    stationaryParticleIBM.C

\*---------------------------------------------------------------------------*/

#ifndef stationaryParticleIBM_H
#define stationaryParticleIBM_H

#include "particleIBM.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class pType> class stationaryParticleIBM;
Ostream& operator<<(Ostream&, const stationaryParticleIBM<pType>&);


template<class pType>
class stationaryParticleIBM
:
    public pType
{
// Private data

    //- Size in bytes of the fields
    static const std::size_t sizeofFields_;

    //- Zero velocity
    vector v_ = Zero;


protected:
// Protected Data

    //- Return the surface velocity given a point
    virtual vector v(const vector& pt) const;

public:

    // Constructors

        //- Construct from a dictionary, used as initial construction of
        //  particle, i.e. copy 0
        stationaryParticleIBM
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const label index
        );

        //- Construct from Istream
        stationaryParticleIBM
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields
        );

        //- Construct as a copy with new center (used with cyclic boundaries)
        stationaryParticleIBM
        (
            const stationaryParticleIBM<pType>& p,
            const label copy,
            const vector& center,
            const vector& theta,
            const vector& v,
            const vector& omega,
            const scalar& age
        );

        //- Construct as a copy
        stationaryParticleIBM(const stationaryParticleIBM<pType>& p);


    //- Destructor
    virtual ~stationaryParticleIBM();


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


    // Public functions

        //- Calcluate needed forcing on fluid
        virtual void forcing
        (
            const surfaceVectorField& Uf,
            const surfaceVectorField& Ufold,
            const surfaceVectorField& Sf,
            volVectorField& F
        ) const;

        //- Calculate and return the drag coefficient
        scalar Cd(const scalar& rhoRef, const vector& UInf) const;


    //- Edit

        //- Move particles in space
        virtual void solve
        (
            const scalar& dt
        );

        //- Update weights, cell lists, and mesh
        virtual void update();


    //- Read/Write functions

        // Write fields
        static void writeFields
        (
            const typename Cloud<stationaryParticleIBM<pType>>& c
        );

        // Read fields
        static void readFields
        (
            typename Cloud<stationaryParticleIBM<pType>>& c
        );


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const stationaryParticleIBM<pType>&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End Foam namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include    "stationaryParticleIBM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

