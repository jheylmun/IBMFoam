/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
\*---------------------------------------------------------------------------*/

#include "slipBC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace velocityBCs
{
    defineTypeNameAndDebug(slipBC, 0);

    addToRunTimeSelectionTable
    (
        velocityBC,
        slipBC,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        velocityBC,
        slipBC,
        copy
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityBCs::slipBC::slipBC
(
    const dictionary& dict,
    const vectorList& mesh,
    const labelList& shellCells,
    const List<labelList>& Is,
    const List<labelList>& Os,
    const List<scalarList>& ws,
    const scalarList& Ws
)
:
    velocityBC(dict, mesh, shellCells, Is, Os, ws, Ws)
{}


Foam::velocityBCs::slipBC::slipBC(const velocityBC& bc)
:
    velocityBC(bc)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityBCs::slipBC::~slipBC()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velocityBCs::slipBC::forcingOnFluid
(
    const scalar& dT,
    const vectorList& U,
    const vectorList& Uold,
    const vectorList& S,
    const particleIBM& p,
    volVectorField& F
) const
{
    forAll(shellCells_, celli)
    {
        label i = shellCells_[celli];

        vector Fi = Zero;
        for(label pti = 0; pti < 4; pti++)
        {
            //  Set inner forcing so that u = -u_outer
            Fi +=
                ws_[celli][pti]/Ws_[celli]
               *(
                    (
                        p.v(mesh_[Os_[celli][pti]])
                      - U[Os_[celli][pti]]
                      - Uold[Os_[celli][pti]]
                    )/dT
                  + S[Os_[celli][pti]]
                );
        }
        for(label pti = 4; pti < 8; pti++)
        {
            //  Set forcing so that u = u_p
            Fi +=
                ws_[celli][pti]/Ws_[celli]
               *(
                    (
                        p.v(mesh_[Is_[celli][pti]])
                      - Uold[Is_[celli][pti]]
                    )/dT
                  + S[Is_[celli][pti]]
                );
        }
        F[i] += Fi;
    }
}
// ************************************************************************* //
