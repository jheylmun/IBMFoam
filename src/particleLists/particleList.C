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

#include "particleList.H"
#include "polyMesh.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"


// * * * * * * * * * * * * Protected Member functions  * * * * * * * * * * * //
template<class pType>
Foam::IOobject Foam::IBM::particleList<pType>::fieldIOobject
(
    const word& fileName,
    const IOobject::readOption& r
) const
{
    return IOobject
    (
        fileName,
        mesh_.time().timePath()/"IBM",
        mesh_,
        r,
        IOobject::AUTO_WRITE,
        false
    );
}

template<class pType>
Foam::tmp<Foam::volSymmTensorField>
Foam::IBM::particleList<pType>::devRhoReff() const
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

template<class pType>
Foam::tmp<Foam::volScalarField> Foam::IBM::particleList<pType>::mu() const
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

template<class pType>
Foam::tmp<Foam::volScalarField> Foam::IBM::particleList<pType>::rho() const
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


template<class pType>
void Foam::IBM::particleList<pType>::computeCollisions()
{
    forAll(*this, i)
    {
        particle& p1 =(*this)[i];
        const vector& x1 = p1.CoM();
        vector& v1 = p1.v();
        scalar mass1 = p1.mass();

        for (label j = i + 1; j < (*this).size(); j++)
        {
            particle& p2 =(*this)[j];
            const vector& x2 = p2.CoM();
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


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class pType>
Foam::IBM::particleList<pType>::particleList
(
    const fvMesh& mesh
)
:
    regIOobject
    (
        IOobject
        (
            "IBMparticles",
            mesh.time().timeName(),
            cloud::prefix,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    ),
    PtrList<pType>(),
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
    PtrList<pType>& pList = *this;
    pList.resize(readLabel(dict_.lookup("nParticles")));

    for (label i = 0; i < pList.size(); i++)
    {
        pList.set
        (
            i,
            new pType
            (
                *this,
                dict_.subDict
                (
                    IOobject::groupName
                    (
                        "particle",
                        Foam::name(i)
                    )
                ),
                Uold_,
                S_
            )
        );
    }
}


template<class pType>
Foam::IBM::particleList<pType>::~particleList()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

template<class pType>
Foam::tmp<Foam::volVectorField> Foam::IBM::particleList<pType>::forcing
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
    PtrList<pType>& pList = *this;

    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    Uold_ = U.oldTime();

    surfaceVectorField Uf(fvc::interpolate(U));
    surfaceVectorField Ufold(fvc::interpolate(U.oldTime()));
    surfaceVectorField Sf(fvc::interpolate(S));

    forAll(pList, particlei)
    {
        pList[particlei].forcing(Uf, Ufold, Sf, F);
    }

    return tmpF;

}

template<class pType>
void Foam::IBM::particleList<pType>::integrateSurfaceStresses()
{
    PtrList<pType>& pList = *this;

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

    forAll(pList, particlei)
    {
        pList[particlei].integrateSurfaceStress(tauf,pf);
    }
}

template<class pType>
void Foam::IBM::particleList<pType>::Cd(IOField<scalar>& Cds) const
{
    const PtrList<pType>& pList = *this;

    forAll(pList, particlei)
    {
        Cds[particlei] = pList[particlei].Cd();
    }
}

template<class pType>
bool Foam::IBM::particleList<pType>::write() const
{
    label np = (*this).size();

    IOField<vector> pos(fieldIOobject("position",IOobject::NO_READ),np);
    IOField<vector> v(fieldIOobject("v",IOobject::NO_READ),np);

    forAll(*this, i)
    {
        const particle& p =(*this)[i];

        pos[i] = p.CoM();
        v[i] = p.v();
    }
    pos.write();
    v.write();

    return true;

}

template<class pType>
bool Foam::IBM::particleList<pType>::writeData(Ostream& os) const
{
    IOField<scalar> cd(fieldIOobject("Cd",IOobject::NO_READ),(*this).size());
    this->Cd(cd);
    forAll(cd, i)
        os << "     Cd." << Foam::name(i) << ": " << cd[i] << endl;
    return os.good();
}


template<class pType>
void Foam::IBM::particleList<pType>::operator++()
{
    computeCollisions();
    integrateSurfaceStresses();
    scalar dt = mesh_.time().deltaT().value();

    vector avgV = Zero;
    scalar avgP = 0.0;
    forAll(*this, i)
    {
        particle& p =(*this)[i];
//         if (!p.onMesh())
//         {
//             p.v() = Zero;
//             continue;
//         }

        p.CoM() += dt*p.v();
        p.v() += dt*p.F()/p.mass();

        p.update();

        avgV += p.v();
        avgP += mag(p.v())*p.mass();
    }

    Info<< "    Average velocity: " << avgV << nl
        << "    Average momentum: " << avgP << endl;
}