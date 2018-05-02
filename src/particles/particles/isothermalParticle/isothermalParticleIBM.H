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
    Foam::isothermalParticleIBM

Description
    abstract class for a particle used with IBM method

SourceFiles
    isothermalParticleIBM.C

\*---------------------------------------------------------------------------*/

#ifndef isothermalParticleIBM_H
#define isothermalParticleIBM_H

#include "particleIBM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

template<class pType> class isothermalParticleIBM;

Ostream& operator<<(Ostream&, const isothermalParticleIBM<pType>&);


template<class pType>
class isothermalParticleIBM
:
    public pType
{
// Private data

    //- Size in bytes of the fields
    static const std::size_t sizeofFields_;


protected:
// Protected Data

    //- Particle Temperature (null)
    scalar T_ = 0.0;


public:

    // Constructors

        //- Construct from a dictionary, used as initial construction of
        //  particle, i.e. copy 0
        isothermalParticleIBM
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const label index
        );

        //- Construct from Istream
        isothermalParticleIBM
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields
        );

        //- Construct as a copy with new center (used with cyclic boundaries)
        isothermalParticleIBM
        (
            const isothermalParticleIBM<pType>& p,
            const label copy,
            const vector& center,
            const vector& theta,
            const vector& v,
            const vector& omega,
            const scalar& age
        );

        //- Construct as a copy
        isothermalParticleIBM(const isothermalParticleIBM<pType>& p);


    //- Destructor
    virtual ~isothermalParticleIBM();


    // Access

        //- Return particle temperature
        virtual const scalar& T() const
        {
            return T_;
        }

        //- Non-const ref to particle Temperature
        virtual scalar& T()
        {
            return T_;
        }


    // Public functions

        //- Calcluate needed forcing on fluid
        virtual void heatFlux
        (
            const surfaceVectorField& Uf,
            const surfaceVectorField& Ufold,
            const surfaceVectorField& Sf,
            volVectorField& F
        ) const;


    //- Edit

//         //- Calculate total stress on particle surface
//         void integrateHeatFlux
//         (
//             const surfaceSymmTensorField& tauf,
//             const surfaceScalarField& p
//         );

        //- Solve particle temperature
        virtual void solve(const scalar& dt);

        //- Update weights, cell lists, and mesh
        virtual void update();


    //- Read/Write functions

        // Write fields
        static void writeFields
        (
            const typename Cloud<isothermalParticleIBM<pType>>& c
        );

        // Read fields
        static void readFields
        (
            typename Cloud<isothermalParticleIBM<pType>>& c
        );


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const isothermalParticleIBM<pType>&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End Foam namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "isothermalParticleIBMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include    "isothermalParticleIBMTmp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

