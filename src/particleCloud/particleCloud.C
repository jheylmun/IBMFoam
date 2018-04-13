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
#include "cyclicFvPatch.H"
#include "processorFvPatch.H"
#include "processorCyclicFvPatch.H"


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

void Foam::particleCloud::combineForces()
{
    forAllIter
    (
        typename Cloud<particleIBM>,
        *this,
        pIter
    )
    {
        pIter().updatedForce() = false;
    }

    if (Pstream::parRun())
    {
        for
        (
            label proci = 0;
            proci < Pstream::nProcs();
            proci++
        )
        {
            Cloud<particleIBM>::const_iterator pIter = (*this).begin();
            boolList end(1, false);
            if (Pstream::myProcNo() == proci)
            {
                if ((*this).cend() == pIter)
                {
                    end = true;
                }
            }
            combineReduce
            (
                end,
                ListPlusEqOp<boolList>()
            );

            while (!end[0])
            {
                labelList pIndex(1, 0);
                boolList pUpdated(1, false);

                if (Pstream::myProcNo() == proci)
                {
                    const particleIBM& p = pIter();
                    pIndex[0] = p.index();
                    pUpdated[0] = p.updatedForce();
                }
                combineReduce
                (
                    pIndex,
                    ListPlusEqOp<labelList>()
                );
                combineReduce
                (
                    pUpdated,
                    ListPlusEqOp<boolList>()
                );

                if (!pUpdated[0])
                {
                    vector totF(Zero);
                    vector totT(Zero);
                    forAllConstIter
                    (
                        typename Cloud<particleIBM>,
                        *this,
                        p2Iter
                    )
                    {
                        if (p2Iter().index() == pIndex[0])
                        {
                            totF += p2Iter().F();
                            totT += p2Iter().T();
                        }
                    }

                    // Combine forces acting on particle pIndex
                    combineReduce
                    (
                        totF,
                        plusEqOp<vector>()
                    );
                    combineReduce
                    (
                        totT,
                        plusEqOp<vector>()
                    );

                    forAllIter
                    (
                        typename Cloud<particleIBM>,
                        *this,
                        p2Iter
                    )
                    {
                        if (p2Iter().index() == pIndex[0])
                        {
                            p2Iter().F() = totF;
                            p2Iter().T() = totT;
                            p2Iter().updatedForce() = true;
                        }
                    }
                }

                end = false;
                if (Pstream::myProcNo() == proci)
                {
                    ++pIter;
                    if ((*this).cend() == pIter)
                    {
                        end = true;
                    }
                }
                combineReduce
                (
                    end,
                    ListPlusEqOp<boolList>()
                );
            }
        }
    }
    else
    {
        forAllConstIter
        (
            typename Cloud<particleIBM>,
            *this,
            pIter
        )
        {
            const particleIBM& p = pIter();
            if (!p.updatedForce())
            {
                vector totF(Zero);
                vector totT(Zero);
                forAllConstIter
                (
                    typename Cloud<particleIBM>,
                    *this,
                    p2Iter
                )
                {
                    if (p2Iter().index() == p.index())
                    {
                        totF += p2Iter().F();
                        totT += p2Iter().T();
                    }
                }

                forAllIter
                (
                    typename Cloud<particleIBM>,
                    *this,
                    p2Iter
                )
                {
                    if (p2Iter().index() == p.index())
                    {
                        p2Iter().F() = totF;
                        p2Iter().T() = totT;
                        p2Iter().updatedForce() = true;
                    }
                }
            }
        }
    }
}


void Foam::particleCloud::setPositions()
{
    if (Pstream::parRun())
    {
        forAllIter
        (
            typename Cloud<particleIBM>,
            *this,
            pIter
        )
        {
            particleIBM& p = pIter();
            if (p.centerOnMesh())
            {
                p.active() = true;
            }
            else
            {
                p.active() = false;
            }
            p.updatedPosition() = false;
        }

        for
        (
            label proci = 0;
            proci < Pstream::nProcs();
            proci++
        )
        {
            Cloud<particleIBM>::const_iterator pIter = (*this).begin();
            boolList end(1, false);
            if (Pstream::myProcNo() == proci)
            {
                if ((*this).cend() == pIter)
                {
                    end = true;
                }
            }
            combineReduce
            (
                end,
                ListPlusEqOp<boolList>()
            );

            while (!end[0])
            {
                labelList pIndex(1, 0);
                labelList pCopy(1, 0);
                boolList pActive(1, false);

                if (Pstream::myProcNo() == proci)
                {
                    const particleIBM& p = pIter();
                    pIndex[0] = p.index();
                    pCopy[0] = p.copy();
                    pActive[0] = p.active();

                }
                combineReduce
                (
                    pIndex,
                    ListPlusEqOp<labelList>()
                );
                combineReduce
                (
                    pCopy,
                    ListPlusEqOp<labelList>()
                );
                combineReduce
                (
                    pActive,
                    ListPlusEqOp<boolList>()
                );

                if (pActive[0])
                {
                    vector pCenter(Zero);
                    vector pTheta(Zero);
                    vector pV(Zero);
                    vector pOmega(Zero);

                    if (Pstream::myProcNo() == proci)
                    {
                        const particleIBM& p = pIter();
                        pCenter = p.center();
                        pTheta = p.theta();
                        pV = p.v();
                        pOmega = p.omega();
                    }

                    combineReduce
                    (
                        pCenter,
                        plusEqOp<vector>()
                    );
                    combineReduce
                    (
                        pTheta,
                        plusEqOp<vector>()
                    );
                    combineReduce
                    (
                        pV,
                        plusEqOp<vector>()
                    );
                    combineReduce
                    (
                        pOmega,
                        plusEqOp<vector>()
                    );

                    forAllIter
                    (
                        typename Cloud<particleIBM>,
                        *this,
                        p2Iter
                    )
                    {
                        if
                        (
                            p2Iter().index() == pIndex[0]
                         && !p2Iter().updatedPosition()
                        )
                        {
                            if (p2Iter().copy() == pCopy[0])
                            {
                                p2Iter().center() = pCenter;
                            }

                            p2Iter().theta() = pTheta;
                            p2Iter().v() = pV;
                            p2Iter().omega() = pOmega;
                        }
                    }
                }
                end = false;
                if (Pstream::myProcNo() == proci)
                {
                    ++pIter;
                    if ((*this).cend() == pIter)
                    {
                        end = true;
                    }
                }
                combineReduce
                (
                    end,
                    ListPlusEqOp<boolList>()
                );
            }
        }
    }
    else
    {
        forAllConstIter
        (
            typename Cloud<particleIBM>,
            *this,
            pIter
        )
        {
            const particleIBM& p = pIter();
            if (!p.active())
            {
                forAllIter
                (
                    typename Cloud<particleIBM>,
                    *this,
                    p2Iter
                )
                {
                    if
                    (
                        p2Iter().index() == p.index()
                     && p2Iter().copy() != p.copy()
                    )
                    {
                        p2Iter().theta() = p.theta();
                        p2Iter().v() = p.v();
                        p2Iter().omega() = p.omega();
                    }
                }
            }
        }
    }
}


void Foam::particleCloud::computeCollisions(const scalar& deltaT)
{
    forAllIter(typename Cloud<particleIBM>, *this, pIter)
    {
        particleIBM& p1 = pIter();

        const vector& x1 = p1.center();
        vector& v1 = p1.v();
        scalar mass1 = p1.mass();

        forAllIter(typename Cloud<particleIBM>, *this, p2Iter)
        {
            particleIBM& p2 = p2Iter();
            const vector& x2 = p2.center();
            vector& v2 = p2.v();
            scalar mass2 = p2.mass();

            vector x12 = x1 - x2;
            scalar magX12 = mag(x12);
            vector v12 = v1 - v2;

            scalar r1 = p1.r(x2);
            scalar r2 = p2.r(x1);

            if (p1.index() > p2.index() && p1.onMesh() && p2.onMesh())
            {
                scalar normalOverlapMag = (r1 + r2) - magX12;

                if (normalOverlapMag > 0 && (x12 & v12) < 0 && (&p1 != &p2))
                {
                    scalar M = mass1*mass2/(mass1 + mass2);

                    vector norm = x12/mag(x12);

                    vector Fn = Zero;
                    vector Ft = Zero;

                    scalar x12Dotv12 = x12 & v12;
                    scalar x12Sqr = magSqr(x12);

                    vector F = 2.0*M*x12Dotv12/x12Sqr*x12*e_/deltaT;
                    Fn = (F & norm)*norm;
                    Ft = F - Fn;

                    p1.Fc()[p2.index()] = -Fn;
                    p2.Fc()[p1.index()] = Fn;
                }
            }
        }
    }

    if (Pstream::parRun())
    {
        List<vectorList> Fc(nParticles_, vectorList(nParticles_, Zero));
        forAllConstIter(typename Cloud<particleIBM>, *this, pIter)
        {
            Fc[pIter().index()] = pIter().Fc();
        }
        combineReduce
        (
            Fc,
            ListListMaxVectorEqOp<List<vectorList>>()
        );

        forAllIter(typename Cloud<particleIBM>, *this, pIter)
        {
            pIter().Fc() = Fc[pIter().index()];
        }
    }
    else
    {
        forAllIter(typename Cloud<particleIBM>, *this, pIter)
        {
            forAllIter(typename Cloud<particleIBM>, *this, p2Iter)
            {
                if (pIter().index() == p2Iter().index())
                {
                    if
                    (
                        mag(pIter().Fc()[p2Iter().index()])
                      > mag(p2Iter().Fc()[pIter().index()])
                    )
                    {
                        p2Iter().Fc()[p2Iter().index()] =
                            pIter().Fc()[pIter().index()];
                    }
                    else
                    {
                        pIter().Fc()[p2Iter().index()] =
                            p2Iter().Fc()[pIter().index()];
                    }
                }
            }
        }
    }
}


bool Foam::particleCloud::found
(
    const label index,
    const label copy,
    const vector& center,
    label& maxCopy
) const
{
    bool f = false;
    forAllConstIter
    (
        typename Cloud<particleIBM>,
        *this,
        pIter
    )
    {
        if
        (
            (
                index == pIter().index()
             && (copy + 1) <= pIter().copy()
            )
         || (
                index == pIter().index()
             && mag(center - pIter().center()) < SMALL
            )
        )
        {
            maxCopy = max(maxCopy, pIter().copy());
            f = true;
        }
    }

    return f;
}


const Foam::particleIBM& Foam::particleCloud::origParticle(const label index) const
{
    forAll(origParticles_, particlei)
    {
        if (index == origParticles_[particlei].index())
        {
            return origParticles_[particlei];
        }
    }
    FatalErrorInFunction
        << "Particle." << index << "not found"
        << exit(FatalError);
    return origParticles_[0];
}


void Foam::particleCloud::computeWallHits(const scalar& deltaT)
{
    if (!walls_.size())
    {
        return;
    }

    forAllIter(typename Cloud<particleIBM>, *this, pIter)
    {
        forAll(walls_, patchi)
        {
            pIter().wallHit(mesh_, deltaT, patchi);
        }
    }

    forAll(walls_, patchi)
    {
        forAllIter(typename Cloud<particleIBM>, *this, pIter)
        {
            forAllIter(typename Cloud<particleIBM>, *this, p2Iter)
            {
                if (pIter().index() == p2Iter().index())
                {
                    if (mag(pIter().Fw()[patchi]) > mag(p2Iter().Fw()[patchi]))
                    {
                        p2Iter().Fw()[patchi] = pIter().Fw()[patchi];
                    }
                    else
                    {
                        pIter().Fw()[patchi] = p2Iter().Fw()[patchi];
                    }
                }
            }
        }
    }
    if (Pstream::parRun())
    {
        List<vectorList> Fw
        (
            nParticles_,
            vectorList(walls_.size(), Zero)
        );
        forAllConstIter(typename Cloud<particleIBM>, *this, pIter)
        {
            Fw[pIter().index()] = pIter().Fw();
        }

        combineReduce
        (
            Fw,
            ListListMaxVectorEqOp<List<vectorList>>()
        );

        forAllIter(typename Cloud<particleIBM>, *this, pIter)
        {
            pIter().Fw() = Fw[pIter().index()];
        }
    }
}


void Foam::particleCloud::computeCyclicHits()
{
    if (!walls_.size())
    {
        return;
    }

    // List of new centers and info on cyclicFvPatch
    // Info = {index, copy, maxCopy}
    List<vectorList> newCyclicVectors; // center, theta, v, omega
    scalarList newCyclicAges;
    labelListList newCyclicInfo;

    // List of centers and info on processorFvPatch
    // Info = {index, copy, neiProc}
    List<vectorList> newProcessorVectors; // center, theta, v, omega
    scalarList newProcessorAges;
    labelListList newProcessorInfo;

    // List of new centers and info on cyclicProcessorFvPatch
    // Info = {index, copy, neiProcID, neiFace, neiProc}
    List<vectorList> newProcessorCyclicVectors; // R, theta, v, omega
    scalarList newProcessorCyclicAges;
    labelListList newProcessorCyclicInfo;

    forAll(patches_, i)
    {
        label patchi = patches_[i];
        forAllConstIter(typename Cloud<particleIBM>, *this, pIter)
        {
            vector R = Zero;
            const particleIBM& p = pIter();
            label facei = p.patchHit
            (
                mesh_,
                patchi,
                R
            );

            if (facei != -1)
            {
                label maxCopy = p.copy();
                if (isA<cyclicFvPatch>(mesh_.boundary()[patchi]))
                {
                    const cyclicPolyPatch& patch =
                        refCast<const cyclicPolyPatch>
                        (
                            mesh_.boundaryMesh()[patchi]
                        );
                    vector newCenter =
                        R
                      + mesh_.Cf().boundaryField()
                        [patch.neighbPatchID()][facei];

                    if (!found(p.index(), p.copy(), newCenter, maxCopy))
                    {
                        newCyclicVectors.append
                        (
                            {
                                newCenter,
                                p.theta(),
                                p.v(),
                                p.omega()
                            }
                        );
                        newCyclicAges.append(p.age());
                        newCyclicInfo.append
                        (
                            labelList({p.index(), p.copy(), maxCopy})
                        );
                    }
                }
                else if (isA<processorFvPatch>(mesh_.boundary()[patchi]))
                {
                    const processorPolyPatch& patch =
                        refCast<const processorPolyPatch>
                        (
                            mesh_.boundaryMesh()[patchi]
                        );

                    newProcessorVectors.append
                    (
                        {
                            p.center(),
                            p.theta(),
                            p.v(),
                            p.omega()
                        }
                    );
                    newProcessorAges.append(p.age());
                    newProcessorInfo.append
                    (
                        labelList({p.index(), p.copy(), patch.neighbProcNo()}));
                }
                else if
                (
                    isA<processorCyclicFvPatch>(mesh_.boundary()[patchi])
                )
                {
                    const processorCyclicPolyPatch& patch =
                        refCast<const processorCyclicPolyPatch>
                        (
                            mesh_.boundaryMesh()[patchi]
                        );
                    newProcessorCyclicVectors.append
                    (
                        {
                            p.center(),
                            p.theta(),
                            p.v(),
                            p.omega()
                        }
                    );
                    newProcessorCyclicAges.append(p.age());
                    newProcessorCyclicInfo.append
                    (
                        labelList
                        (
                            {
                                p.index(),
                                p.copy(),
                                patch.referPatchID(),
                                facei,
                                patch.neighbProcNo()
                            }
                        )
                    );
                }
            }
        }
    }

    forAll(newCyclicInfo, newParticlei)
    {
        this->append
        (
            new particleIBM
            (
                origParticle(newCyclicInfo[newParticlei][0]),
                newCyclicInfo[newParticlei][2] + 1,
                newCyclicVectors[newParticlei][0],
                newCyclicVectors[newParticlei][1],
                newCyclicVectors[newParticlei][2],
                newCyclicVectors[newParticlei][3],
                newCyclicAges[newParticlei]
            )
        );
    }

    if (Pstream::parRun())
    {
        combineReduce
        (
            newProcessorVectors,
            combineListEqOp<List<vectorList>>()
        );
        combineReduce
        (
            newProcessorAges,
            combineListEqOp<scalarList>()
        );
        combineReduce
        (
            newProcessorInfo,
            combineListEqOp<labelListList>()
        );
        combineReduce
        (
            newProcessorCyclicVectors,
            combineListEqOp<List<vectorList>>()
        );
        combineReduce
        (
            newProcessorCyclicAges,
            combineListEqOp<scalarList>()
        );
        combineReduce
        (
            newProcessorCyclicInfo,
            combineListEqOp<labelListList>()
        );

        forAll(newProcessorInfo, newParticlei)
        {
            if
            (
                Pstream::myProcNo()
             == newProcessorInfo[newParticlei][2]
            )
            {
                label index = newProcessorInfo[newParticlei][0];
                label copy = newProcessorInfo[newParticlei][1];
                vector center = newProcessorVectors[newParticlei][0];
                label maxCopy = 0;
                if (!found(index, copy, center, maxCopy))
                {
                    this->append
                    (
                        new particleIBM
                        (
                            origParticle(index),
                            maxCopy + 1,
                            center,
                            newProcessorVectors[newParticlei][1],
                            newProcessorVectors[newParticlei][2],
                            newProcessorVectors[newParticlei][3],
                            newProcessorAges[newParticlei]
                        )
                    );
                }
            }
        }

        forAll(newProcessorCyclicInfo, newParticlei)
        {
            if
            (
                Pstream::myProcNo() == newProcessorCyclicInfo[newParticlei][4]
            )
            {
                label index = newProcessorCyclicInfo[newParticlei][0];
                label copy = newProcessorCyclicInfo[newParticlei][1];
                vector newCenter =
                    newProcessorCyclicVectors[newParticlei][0]
                  + mesh_.Cf().boundaryField()
                    [newProcessorCyclicInfo[newParticlei][2]][newProcessorCyclicInfo[newParticlei][3]];

                label maxCopy = 0;
                if (!found(index, copy, newCenter, maxCopy))
                {
                    this->append
                    (
                        new particleIBM
                        (
                            origParticle(index),
                            maxCopy + 1,
                            newCenter,
                            newProcessorCyclicVectors[newParticlei][1],
                            newProcessorCyclicVectors[newParticlei][2],
                            newProcessorCyclicVectors[newParticlei][3],
                            newProcessorCyclicAges[newParticlei]
                        )
                    );
                }
            }
        }
    }
}


void Foam::particleCloud::deleteInActiveParticles()
{
    forAllIter(typename Cloud<particleIBM>, *this, pIter)
    {
        particleIBM& p = pIter();
        if (!p.onMesh())
        {
            this->deleteParticle(p);
        }
    }
}

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

Foam::particleCloud::particleCloud
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    Cloud<particleIBM>
    (
        mesh,
        "IBMparticles",
        false
    ),
    mesh_(mesh),
    dict_(dict),
    nParticles_(readLabel(dict_.lookup("nParticles"))),
    origParticles_(nParticles_),
    UInf_
    (
        dict_.subDict("flow").template lookupOrDefault<vector>("UInf", Zero)
    ),
    flowNormal_(UInf_/(mag(UInf_) + SMALL)),
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
    coupled_(dict_.lookup("coupled")),
    moving_(dict_.lookup("moving")),
    rotating_(dict_.lookup("rotating"))
{
    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            isA<cyclicFvPatch>(mesh_.boundary()[patchi])
         || isA<processorFvPatch>(mesh_.boundary()[patchi])
         || isA<processorCyclicFvPatch>(mesh_.boundary()[patchi])
        )
        {
            patches_.append(patchi);
        }
        else if (isA<wallFvPatch>(mesh_.boundary()[patchi]))
        {
            walls_.append(patchi);
        }
    }

    for (label i = 0; i < nParticles_; i++)
    {
        const dictionary& pDict =
            dict_.subDict
            (
                IOobject::groupName
                (
                    "particle",
                    Foam::name(i)
                )
            );

        this->append
        (
            new particleIBM
            (
                mesh_,
                pDict,
                i
            )
        );
    }

    label i = 0;
    forAllConstIter(typename Cloud<particleIBM>, *this, pIter)
    {
        const particleIBM& p = pIter();
        origParticles_.set(i, new particleIBM(p));
        origParticles_[i].setCollisions(nParticles_, walls_.size());
        i++;
    }

    forAllIter(typename Cloud<particleIBM>, *this, pIter)
    {
        pIter().setCollisions(nParticles_, walls_.size());
    }

    deleteInActiveParticles();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

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

    if (!coupled_)
    {
        return tmpF;
    }

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
        i++;
    }
}

void Foam::particleCloud::solve()
{
    integrateSurfaceStresses();
    combineForces();

    scalar dt = mesh_.time().deltaT().value();
    scalar localDt = dt/10.0;
    scalar localTime = localDt;

    while (localTime < dt)
    {
        if (moving_)
        {
            computeCollisions(localDt);
            computeWallHits(localDt);
            computeCyclicHits();
        }
        combineForces();

        forAllIter(typename Cloud<particleIBM>, *this, pIter)
        {
            particleIBM& p = pIter();
            p.solve(localDt, moving_, rotating_);
        }
        setPositions();

        if (localTime + localDt > dt)
        {
            localDt = dt - localTime;
        }
        else
        {
            localTime += localDt;
        }
    }

    vector avgV = Zero;
    vector avgOmega = Zero;
    vector avgP = Zero;
    scalarList totK(1, 0.0);
    forAllConstIter(typename Cloud<particleIBM>, *this, pIter)
    {
        const particleIBM& p = pIter();

        avgOmega += p.omega()*p.active();
        avgV += p.v()*p.active();
        avgP += p.v()*p.mass()*p.active();
        totK[0] += 0.5*p.mass()*(p.v() & p.v())*p.active();
    }

    forAllIter
    (
        typename Cloud<particleIBM>,
        *this,
        pIter
    )
    {
        pIter().update();
    }
    deleteInActiveParticles();

    label nParticles = Cloud<particleIBM>::size();
    labelList totParticles(1, nParticles);
    combineReduce
    (
        totParticles,
        ListPlusEqOp<labelList>()
    );

    if (nParticles)
    {
        combineReduce(avgV, plusEqOp<vector>());
        combineReduce(avgOmega, plusEqOp<vector>());
        combineReduce(avgP, plusEqOp<vector>());
        combineReduce(totK, ListPlusEqOp<scalarList>());

        Info<< "    Actual number of active particles: " << nParticles_ << nl
            << "    Total number of active particles: " << totParticles << nl
            << "    Average velocity: " << avgV/nParticles << nl
            << "    Average omega: " << avgOmega/nParticles << nl
            << "    Average momentum: " << avgP/nParticles << nl
            << "    total K: " << totK[0] << endl;
    }

    IOField<scalar> cds(fieldIOobject("Cd", IOobject::NO_READ), this->size());
    Cd(cds);
}

// void Foam::particleCloud::computeCyclicHits()
// {
//     if (!walls_)
//     {
//         return;
//     }
//
//     vectorList newCenters;
//     labelListList newParticleInfo;
//
//     forAll(mesh_.boundary(), patchi)
//     {
//         if (isA<cyclicFvPatch>(mesh_.boundary()[patchi]))
//         {
//             const cyclicPolyPatch& patch =
//                 refCast<const cyclicPolyPatch>
//                 (
//                     mesh_.boundaryMesh()[patchi]
//                 );
//             label neiPatchID = patch.neighbPatchID();
//
//             forAllConstIter(typename Cloud<particleIBM>, *this, pIter)
//             {
//                 const particleIBM& p = pIter();
//                 vector R = Zero;
//                 label facei = p.cyclicHit
//                 (
//                     mesh_,
//                     patchi,
//                     R
//                 );
//
//                 if (facei != -1)
//                 {
//                     vector newCenter =
//                         R + mesh_.Cf().boundaryField()[neiPatchID][facei];
//
//                     bool found = false;
//                     label maxCopy = p.copy();
//                     forAllConstIter
//                     (
//                         typename Cloud<particleIBM>,
//                         *this,
//                         p2Iter
//                     )
//                     {
//                         if
//                         (
//                             p.index() == p2Iter().index()
//                          && (maxCopy + 1) <= p2Iter().copy()
//                         )
//                         {
//                             found = true;
//                             maxCopy = max(maxCopy, p2Iter().copy());
//                         }
//                     }
//
//                     if (!found)
//                     {
//                         newCenters.append(newCenter);
//                         newParticleInfo.append
//                         (
//                             labelList({p.index(), p.copy(), maxCopy})
//                         );
//                     }
//                 }
//             }
//         }
//     }
//
//     forAll(newCenters, newParticlei)
//     {
//         forAllConstIter
//         (
//             typename Cloud<particleIBM>,
//             *this,
//             pIter
//         )
//         {
//             if
//             (
//                 pIter().index() == newParticleInfo[newParticlei][0]
//              && pIter().copy() == newParticleInfo[newParticlei][1]
//             )
//             {
//                 this->append
//                 (
//                     new particleIBM
//                     (
//                         pIter(),
//                         newCenters[newParticlei],
//                         newParticleInfo[newParticlei][2] + 1
//                     )
//                 );
//             }
//         }
//     }
// }