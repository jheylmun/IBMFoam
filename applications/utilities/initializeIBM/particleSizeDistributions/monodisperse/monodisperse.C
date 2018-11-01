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
    ds_(pTypes_.size(), 0.0)
{
    scalar totalV = 0.0;
    forAll(mesh_.cellVolumes(), celli)
    {
        totalV += mesh_.cellVolumes()[celli];
    }

    label nDims = mesh_.nGeometricD();
    scalar totParticalVolumeFraction(0.0);

    forAll(pTypes_, pTypei)
    {
        const dictionary& pDict = dict_.subDict(pTypes_[pTypei]);
        ds_[pTypei] = readScalar(pDict.lookup("d"));

        scalar volumeFractioni = pDict.lookupOrDefault("volumeFraction", 0.0);
        nParticles_[pTypei] = pDict.lookupOrDefault("nParticles", 0);

        scalar oneParticleV(0.0);
        if (nDims == 3)
        {
            oneParticleV = Foam::constant::mathematical::pi/6.0*pow3(ds_[pTypei]);
        }
        else
        {
            oneParticleV =
                Foam::constant::mathematical::pi/4*sqr(ds_[pTypei])
            *(bb_.max().z() - bb_.min().z());
        }

        if (nParticles_[pTypei] == 0 && volumeFractioni == 0.0)
        {
            FatalErrorInFunction
                << "Either the number of particles or the volume fraction " << nl
                << " must be specified to initialize a field of IBM particles."
                << nl << exit(FatalError);
        }

        if (volumeFractioni > 0)
        {
            nParticles_[pTypei] = label(volumeFractioni*totalV/oneParticleV);
        }
        else
        {
            volumeFractioni = nParticles_[pTypei]*oneParticleV/totalV;
        }
        nTotParticles_ += nParticles_[pTypei];
        totParticalVolumeFraction += volumeFractioni;

        for  (label i = 0; i < nParticles_[pTypei]; i++)
        {
            pShapes_.append
            (
                particleShape::New(mesh_, pDict, Zero, false).ptr()
            );
        }
    }

    Info<< "Created " << nTotParticles_ << " particles." << nl
        << "Average volume fraction = " << totParticalVolumeFraction << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IBM::monodisperse::~monodisperse()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IBM::monodisperse::writeParticles(Ostream& os)
{
    label threeD = (mesh_.nGeometricD() > 2 ? 1 : 0);

    label particlei = 0;
    forAll(pTypes_, pTypei)
    {
        Info<< "creating " << nParticles_[pTypei] << " "
            << pTypes_[pTypei] << " particles" << endl;

        const dictionary& pDict = dict_.subDict(pTypes_[pTypei]);

        for (label pi = 0; pi < nParticles_[pTypei]; pi++)
        {
            os  << IOobject::groupName("particle", Foam::name(particlei)) << nl
                << token::BEGIN_BLOCK <<endl;

            vector pos = Zero;

            bool validPos = false;
            while (!validPos)
            {
                validPos = true;
                vector diff(bb_.max() - bb_.min());
                pos =
                    bb_.min()
                  + vector
                    (
                        rand_.scalar01()*diff.x(),
                        rand_.scalar01()*diff.y(),
                        rand_.scalar01()*diff.z()
                    );
                if (!threeD)
                {
                    pos.z() = 0.5*(bb_.min().z() + bb_.max().z());
                }
                pShapes_[particlei].center() = pos;

                for (label i = 0; i < particlei; i++)
                {
                    if
                    (
                        (
                            mag(pos - pShapes_[i].center())
                         <= (
                                pShapes_[particlei].r(pShapes_[i].center())
                              + pShapes_[i].r(pos)
                            )
                        )
                     || this->mesh_.findCell(pos) == -1
                    )
                    {
                        validPos = false;
                        break;
                    }
                }
            }

            os  << "particleShape" << token::TAB
                << pDict.lookupType<word>("particleShape")
                << token::END_STATEMENT << nl

                << "position" << token::TAB << pShapes_[particlei].center()
                << token::END_STATEMENT << nl

                << pShapes_[particlei]

                << "nTheta" << token::TAB
                << readLabel(pDict.lookup("nTheta")) << token::END_STATEMENT << nl

                << "nk" << token::TAB
                << pDict.lookupOrDefault("nk", 1) << token::END_STATEMENT << nl

                << "delta" << token::TAB
                << readScalar(pDict.lookup("delta")) << token::END_STATEMENT << nl;

            Switch randRho = pDict.lookupOrDefault("randomRho", false);
            scalar rhoInit = pDict.lookupOrDefault("rho", 0.0);
            if (randRho)
            {
                autoPtr<randomDistribution> dist
                (
                    randomDistribution::New
                    (
                        pDict.subDict("rhoDistributionCoeffs")
                    )
                );

                rhoInit = dist->RV(rand_);
            }
            os  << "rho" << token::TAB
                << rhoInit << token::END_STATEMENT << nl;

                Switch randTheta = pDict.lookupOrDefault("randomTheta", false);
            vector thetaInit = pDict.lookupOrDefault<vector>("theta", Zero);
            if (randTheta)
            {
                autoPtr<randomDistribution> dist
                (
                    randomDistribution::New
                    (
                        pDict.subDict("thetaDistributionCoeffs")
                    )
                );

                thetaInit.z() = dist->RV(rand_);
                if (threeD)
                {
                    thetaInit.x() = dist->RV(rand_);
                    thetaInit.y() = dist->RV(rand_);
                }
            }
            if (mag(thetaInit) > 0)
            {
                os  << "theta" << token::TAB
                    << thetaInit << token::END_STATEMENT << nl;
            }

            Switch randV = pDict.lookupOrDefault("randomVelocity", false);
            vector vInit = pDict.lookupOrDefault<vector>("v", Zero);
            if (randV)
            {
                autoPtr<randomDistribution> dist
                (
                    randomDistribution::New
                    (
                        pDict.subDict("velocityDistributionCoeffs")
                    )
                );

                vInit.x() = dist->RV(rand_);
                vInit.y() = dist->RV(rand_);

                if (threeD)
                {
                    vInit.z() = dist->RV(rand_);
                }
            }
            if (mag(vInit) > 0)
            {
                os  << "v" << token::TAB
                    << vInit << token::END_STATEMENT << nl;
            }

            Switch randOmega = pDict.lookupOrDefault("randomOmega", false);
            vector omegaInit = pDict.lookupOrDefault<vector>("omega", Zero);
            if (randOmega)
            {
                autoPtr<randomDistribution> dist
                (
                    randomDistribution::New
                    (
                        pDict.subDict("omegaDistributionCoeffs")
                    )
                );

                omegaInit.z() = dist->RV(rand_);
                if (threeD)
                {
                    omegaInit.x() = dist->RV(rand_);
                    omegaInit.y() = dist->RV(rand_);
                }
            }
            if (mag(omegaInit) > 0)
            {
                os  << "omega" << token::TAB
                    << omegaInit << token::END_STATEMENT << nl;
            }

            os << token::END_BLOCK <<endl;

            Info<< "Particle " << particlei << " placed." << endl;
            particlei++;
        }
    }
    return os.good();
}
