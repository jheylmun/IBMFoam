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

#include "particleCloud.H"
#include "polyMesh.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "wallFvPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(particleCloud, 0);
}

// * * * * * * * * * * * * Protected Member functions  * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField>
Foam::particleCloud::devRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            mesh_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff();
    }
    else if (mesh_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const incompressible::turbulenceModel& turb =
            mesh_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return rho()*turb.devReff();
    }
    else if (mesh_.foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const fluidThermo& thermo =
            mesh_.lookupObject<fluidThermo>(fluidThermo::dictName);

        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

        return -thermo.mu()*dev(twoSymm(fvc::grad(U)));
    }
    else if
    (
        mesh_.foundObject<transportModel>("transportProperties")
    )
    {
        const transportModel& laminarT =
            mesh_.lookupObject<transportModel>("transportProperties");

        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

        return -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (mesh_.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             mesh_.lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu(transportProperties.lookup("nu"));

        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

        return -rho()*nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}

Foam::tmp<Foam::volScalarField> Foam::particleCloud::mu() const
{
    if (mesh_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo =
             mesh_.lookupObject<fluidThermo>(basicThermo::dictName);

        return thermo.mu();
    }
    else if
    (
        mesh_.foundObject<transportModel>("transportProperties")
    )
    {
        const transportModel& laminarT =
            mesh_.lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (mesh_.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             mesh_.lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu
        (
            "nu",
            dimViscosity,
            transportProperties.lookup("nu")
        );

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}

Foam::tmp<Foam::volScalarField> Foam::particleCloud::rho() const
{
    if (!mesh_.foundObject<volScalarField>("thermo:rho"))
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("rho", dimDensity, rhoRef_)
            )
        );
    }
    else
    {
        return(mesh_.lookupObject<volScalarField>("thermo:rho"));
    }
}


void Foam::particleCloud::computeCollisions()
{
    forAllIter(typename Cloud<particleIBM>, *this, pIter)
    {
        particleIBM& p1 = pIter();

        const vector& x1 = p1.center();
        vector& v1 = p1.v();
        scalar mass1 = p1.mass();

        forAllIter(typename Cloud<particleIBM>, *this, pIter)
        {
            particleIBM& p2 = pIter();
            const vector& x2 = p2.center();
            vector& v2 = p2.v();
            scalar mass2 = p2.mass();

            vector x12 = x1 - x2;
            vector v12 = v1 - v2;

            scalar r1 = p1.r(x2);
            scalar r2 = p2.r(x1);

            if (mag(x12) <= (r1 + r2) && (x12 & v12) < 0)
            {
                scalar M = mass1 + mass2;
                scalar x12Dotv12 = x12 & v12;
                scalar x12Sqr = magSqr(x12);

                v1 -= 2.0*mass2/M*x12Dotv12/x12Sqr*x12*e_;
                v2 += 2.0*mass1/M*x12Dotv12/x12Sqr*x12*e_;
            }
        }
    }
}


void Foam::particleCloud::computeWallHits()
{
    if (walls_)
    {
        return;
    }

    forAllIter(typename Cloud<particleIBM>, *this, pIter)
    {
        particleIBM& p = pIter();
        p.wallHit(mesh_);
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::particleCloud::particleCloud
(
    const fvMesh& mesh
)
:
    Cloud<particleIBM>
    (
        mesh,
        "IBMparticles",
        false
    ),
    mesh_(mesh),
    dict_
    (
        IOobject
        (
            "initialStates",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    UInf_
    (
        dict_.subDict("flow").template lookupOrDefault<vector>("UInf", Zero)
    ),
    flowNormal_(UInf_/(mag(UInf_) + small)),
    rhoRef_
    (
        dict_.subDict("flow").template lookupOrDefault<scalar>("rhoInf", 1.0)
    ),
    pRef_
    (
        dict_.subDict("flow").template lookupOrDefault<scalar>("pInf", 0.0)
    ),
    e_(readScalar(dict_.lookup("e"))),
    Uold_
    (
        IOobject
        (
            "Uold",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, Zero)
    ),
    S_
    (
        IOobject
        (
            "ibmSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        mesh,
        dimensionedVector("zero", dimAcceleration, Zero)
    ),
    moving_(dict_.lookup("moving"))
{
    label nParticles = readLabel(dict_.lookup("nParticles"));
    for (label i = 0; i < nParticles; i++)
    {
        this->append
        (
            new particleIBM
            (
                mesh_,
                dict_.subDict
                (
                    IOobject::groupName
                    (
                        "particle",
                        Foam::name(i)
                    )
                )
            )
        );
    }

    forAll(mesh_.boundary(), patchi)
    {
        if (isA<wallFvPatch>(mesh_.boundary()[patchi]))
        {
            walls_ = true;
            continue;
        }
    }
}


Foam::particleCloud::~particleCloud()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::particleCloud::forcing
(
    const volVectorField& S
)
{
    tmp<volVectorField> tmpF
    (
        new volVectorField
        (
            IOobject
            (
                "F",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector
            (
                "F",
                dimForce/dimVolume/dimDensity,
                Zero
            )
        )
    );

    volVectorField& F = tmpF.ref();

    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    Uold_ = U.oldTime();

    surfaceVectorField Uf(fvc::interpolate(U));
    surfaceVectorField Ufold(fvc::interpolate(U.oldTime()));
    surfaceVectorField Sf(fvc::interpolate(S));

    forAllConstIter(typename Cloud<particleIBM>, *this, pIter)
    {
        const particleIBM& p = pIter();
        p.forcing(Uf, Ufold, Sf, F);
    }

    return tmpF;

}

void Foam::particleCloud::integrateSurfaceStresses()
{
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    dimensionedScalar rho("rho", dimDensity, rhoRef());
    dimensionedScalar pRef =
        dimensionedScalar("pRef", dimPressure, pRef_)/rho;
    if (p.dimensions() != dimPressure)
    {
        pRef.dimensions().reset(dimPressure/dimDensity);
    }

    surfaceSymmTensorField tauf(fvc::interpolate(devRhoReff()));
    surfaceScalarField pf(fvc::interpolate(p) - pRef);

    forAllIter(typename Cloud<particleIBM>, *this, pIter)
    {
        particleIBM& p = pIter();
        p.integrateSurfaceStress(tauf,pf);
    }
}

void Foam::particleCloud::Cd(IOField<scalar>& Cds) const
{
    label i = 0;
    forAllConstIter(typename Cloud<particleIBM>, *this, pIter)
    {
        const particleIBM& p = pIter();
        Cds[i] = p.Cd(rhoRef_, UInf_);
    }
}

void Foam::particleCloud::operator++()
{
    integrateSurfaceStresses();
    if (!moving_)
    {
        return;
    }
    computeCollisions();
    computeWallHits();
    scalar dt = mesh_.time().deltaT().value();

    vector avgV = Zero;
    scalar avgP = 0.0;
    forAllIter(typename Cloud<particleIBM>, *this, pIter)
    {
        particleIBM& p = pIter();
        p.solve(dt);

        avgV += p.v();
        avgP += mag(p.v())*p.mass();
    }

    Info<< "    Average velocity: " << avgV << nl
        << "    Average momentum: " << avgP << endl;
}