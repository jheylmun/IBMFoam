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

#include "particleSizeDistribution.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace IBM
{
    defineTypeNameAndDebug(particleSizeDistribution, 0);
    defineRunTimeSelectionTable(particleSizeDistribution, dictionary);
}
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IBM::particleSizeDistribution::particleSizeDistribution
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    IOobject
    (
        "PSD",
        mesh.time().constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_(mesh),
    rand_(mesh_.time().elapsedClockTime()),
    dict_(dict),
    nParticles_(dict.subDict("PSD").lookupOrDefault<label>("nParticles",1)),
    volumeFraction_
    (
        dict.subDict("PSD").lookupOrDefault<scalar>
        (
            "volumeFraction",
            0
        )
    ),
    deltaR_(readScalar(dict.subDict("PSD").lookup("deltaR"))),
    nTheta_
    (
        readLabel(dict.subDict("baseParticle").lookup("nTheta"))
    ),
    nk_(dict.subDict("baseParticle").lookupOrDefault<label>("nk",3)),
    delta_(readScalar(dict.subDict("baseParticle").lookup("delta")))
{
    if (nParticles_ == -1 && volumeFraction_ == -1.0)
    {
        FatalErrorInFunction
            << "Either the number of particles or the volume fraction " << nl
            << " must be specified to initialize a field of IBM particles."
            << nl << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IBM::particleSizeDistribution::~particleSizeDistribution()
{}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

bool Foam::IBM::particleSizeDistribution::write()
{
    OFstream os
    (
        mesh_.time().constant()/"initialStates"
    );
    IOdictionary initStateHeader
    (
        IOobject
        (
            "initialStates",
            mesh_.time().constant(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    initStateHeader.writeHeader(os);

    os  << "nParticles" << token::TAB << nParticles_
        << token::END_STATEMENT << endl;

    dictionary flowDict(dict_.subDict("flow"));

    os  << "flow" << nl
        << token::BEGIN_BLOCK << nl
        << "UInf" << token::TAB << vector(flowDict.lookup("UInf"))
        << token::END_STATEMENT << nl;

    if (flowDict.found("rho"))
    {
        os  << "rho" << token::TAB << word(flowDict.lookup("rho"))
            << token::END_STATEMENT << nl
            << "rhoInf" << token::TAB << readScalar(flowDict.lookup("rhoInf"))
            << token::END_STATEMENT << nl;
    }

    if (flowDict.found("pRef"))
    {
        os  << "p" << token::TAB << word(flowDict.lookup("p"))
            << token::END_STATEMENT << nl
            << "pRef" << token::TAB << readScalar(flowDict.lookup("pRef"))
            << token::END_STATEMENT << nl;
    }
    os  << token::END_BLOCK << endl;

    writeParticles(os);

    os << endl;

    return os.good();
}
