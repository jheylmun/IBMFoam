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
#include "Time.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

Foam::IOobject Foam::IBM::particleList::fieldIOobject
(
    const word& fieldName,
    const IOobject::readOption r
) const
{
    IOobject a
    (
        fieldName,
        time().timeName(),
        *this,
        r,
        IOobject::NO_WRITE,
        false
    );
    return a;
}


void Foam::IBM::particleList::writeFields() const
{
    particle::writeFields(*this);
//     const particleList& pList = (*this);
//     label np = pList.size();
//
//     IOField<vector> pos(fieldIOobject("position",IOobject::NO_READ),np);
//     IOField<vector> v(fieldIOobject("v",IOobject::NO_READ),np);
//
//     forAll(pList, i)
//     {
//         const particle& p =(*this)[i];
//
//         pos[i] = p.CoM();
//         v[i] = p.v();
//     }
//     pos.write();
//     v.write();
//
//     return true;

}

bool Foam::IBM::particleList::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    const PtrList<particle>& pList = *this;
    writeFields();
    return cloud::writeObject(fmt, ver, cmp, pList.size());

//     IOField<scalar> cd(fieldIOobject("Cd",IOobject::NO_READ),(*this).size());
//     this->Cd(cd);
//     forAll(cd, i)
//         os << "     Cd." << Foam::name(i) << ": " << cd[i] << endl;
//     return os.good();
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

Foam::Ostream& Foam::IBM::operator<<(Ostream& os, const particleList& pc)
{
    pc.writeData(os);

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const particleList&)");

    return os;
}


// ************************************************************************* //