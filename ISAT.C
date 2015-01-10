/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "ISAT.H"
#include "chemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::ISAT<ChemistryModel>::ISAT
(
    const fvMesh& mesh
)
:
    chemistrySolver<ChemistryModel>(mesh),
    W_(this->nSpecie())
{
//      ncv     - number of composition variables (integer)
//      nfullv  - number of items in the full representation (integer)
//      nstrms  - number of streams (integer)
    int ncv, nfull, nstrms;
    ciinit_( &ncv, &nfull, &nstrms );
    label nSpecie = this->nSpecie();
    if (ncv != nSpecie+1)
    {
        FatalErrorIn("ChemistryModel::ISAT::ISAT(const mesh&)")
            << "The number of species in OpenFOAM does not match "
            << "the number of species in ISAT library:" << nl
            << "Check chem.bin file for errors" << nl
            << abort(FatalError);
    }

    for (register int i=0; i<nSpecie; i++)
    {
        W_[i] = this->specieThermo()[i].W();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::ISAT<ChemistryModel>::~ISAT()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::ISAT<ChemistryModel>::solve
(
    scalarField& c, // kmol/m3
    scalar& T,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{ // W: kg/kmol
//
// Note: Preferrably, we should only call cirxn_ here, but we need
// sensible enthalpy, hs, when calling cirxn_. To get hs, we call ciconv_
// This gives some overhead, which we would like to get rid of... do not know how yet
    int nSpecie = this->nSpecie();
    int ncv = nSpecie+1;
    double Z0[ncv];  // initial specific mole number [mole/g] array + hs
    double Z1[ncv];  // specific mole number after deltaT array + hs
    double X0[nSpecie]; // initial molefrac

    scalar moleTot = 0.;
    double densitySI = 0.;
    for (register int i=0; i<nSpecie; i++)
    {
        densitySI += c[i]*W_[i];
        moleTot += c[i];
        X0[i] = c[i];
    }
    for (register int i=0; i<nSpecie; i++)
    {
        X0[i] /= moleTot;
    }
    double densityCGS=1.;
    double pCGS;
    double hs=-1.;
    int jd=2;
    int jp=3;
    int je=1;
    int js=1;
    int kd=2;
    int kp=2;
    int ke=2;
    int ks=3;
    ciconv_(&ncv, &jd, &jp, &je, &js, &densitySI, &p, &T, X0, &kd, &kp, &ke, &ks, &densityCGS, &pCGS, &hs, Z0);
    Z0[nSpecie] = hs;
    double dpt[3] = {densityCGS, pCGS, -1};
    cirxn_(&deltaT, &ncv, Z0, Z1, dpt);

    for (register int i=0; i<nSpecie; i++)
    {
        c[i] = Z1[i]*densitySI;
    }
    subDeltaT = deltaT;
}


// ************************************************************************* //
