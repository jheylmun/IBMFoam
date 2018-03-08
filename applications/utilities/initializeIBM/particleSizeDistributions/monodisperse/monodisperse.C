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

#include "monodisperse.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace IBM
{
    defineTypeNameAndDebug(monodisperse, 0);

    addToRunTimeSelectionTable
    (
        particleSizeDistribution,
        monodisperse,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IBM::monodisperse::monodisperse
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    particleSizeDistribution(mesh,dict),
    d_(readScalar(dict.subDict("baseParticle").lookup("d"))),
    bb1_(dict.subDict("PSD").lookup("bb1")),
    bb2_(dict.subDict("PSD").lookup("bb2")),
    p1_(dict.subDict("baseParticle").lookup("p1")),
    p2_(dict.subDict("baseParticle").lookup("p2"))
{
    scalar totalV = 0.0;
    forAll(mesh_.cellVolumes(), celli)
    {
        totalV += mesh_.cellVolumes()[celli];
    }

    bool threeD = 1;
    if (mesh_.nGeometricD() == 2) threeD = 0;
    scalar oneParticleV(0.0);

    if (threeD)
    {
        oneParticleV = Foam::constant::mathematical::pi/6.0*pow3(d_);
    }
    else
    {
        oneParticleV =
            Foam::constant::mathematical::pi/4*sqr(d_)
           *mag(p2_ - p1_);
    }

    if (volumeFraction_ != 0)
    {
        nParticles_ = label(volumeFraction_*totalV/oneParticleV);
    }
    else
    {
        volumeFraction_ = nParticles_*oneParticleV/totalV;
    }

    particlePos_ = vectorList(nParticles_, Zero);

    Info<< "Created " << nParticles_ << " particles." << nl
        << "Volume fraction = " << volumeFraction_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IBM::monodisperse::~monodisperse()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IBM::monodisperse::writeParticles(Ostream& os)
{
    bool threeD = 1;
    if (mesh_.nGeometricD() == 2)
    {
        threeD = 0;
    }

    forAll(particlePos_, particlei)
    {
        os  << IOobject::groupName("particle",Foam::name(particlei)) << nl
            << token::BEGIN_BLOCK <<endl;

        vector pos = Zero;

        bool validPos = false;

        while (!validPos)
        {
            validPos = true;
            vector diff(bb2_ - bb1_);
            pos =
                bb1_
              + vector
                (
                    rand_.scalar01()*diff.x(),
                    rand_.scalar01()*diff.y(),
                    rand_.scalar01()*diff.z()*threeD
                );
            for (label i = 0; i < particlei; i++)
            {
                if (mag(pos - particlePos_[i]) - d_ <= deltaR_)
                {
                    validPos = false;
                    break;
                }
            }
        }

        particlePos_[particlei] = pos;

        if (threeD)
        {
            os  << "particleShape" << token::TAB
                << "sphere" << token::END_STATEMENT << nl

                << "center" << token::TAB
                << pos << token::END_STATEMENT << nl

                << "d" << token::TAB
                << d_ << token::END_STATEMENT << nl

                << "nTheta" << token::TAB
                << nTheta_ << token::END_STATEMENT << nl

                << "nk" << token::TAB
                << nk_<< token::END_STATEMENT << nl

                << "delta" << token::TAB
                << delta_ << token::END_STATEMENT << nl;
        }
        else
        {
            os  << "particleShape" << token::TAB
                << "cylinder" << token::END_STATEMENT << nl

                << "p1" << token::TAB << vector(pos.x(),pos.y(),p1_.z())
                << token::END_STATEMENT  << nl

                << "p2" << token::TAB << vector(pos.x(),pos.y(),p2_.z())
                << token::END_STATEMENT << nl

                << "d" << token::TAB
                << d_ << token::END_STATEMENT  << nl

                << "nTheta" << token::TAB
                << nTheta_ << token::END_STATEMENT << nl

                << "delta" << token::TAB
                << delta_ << token::END_STATEMENT << nl;
        }
        os << token::END_BLOCK <<endl;
    }
    return os.good();
}
