/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2014 Karl-Johan Nogenmyr
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
#include "specie.H"
#include "chemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::ISAT<ChemistryModel>::ISAT
(
    const fvMesh& mesh
)
:
    chemistrySolver<ChemistryModel>(mesh),
    coeffsDict_(this->subDict("ISATCoeffs")),
    W_(this->nSpecie()),
    href_(this->nSpecie()),
    ncv_(-1)
{
//      nfullv  - number of items in the full representation (integer)
//      nstrms  - number of streams (integer)
    int nfull, nstrms;

    label nSpecie = this->nSpecie();

    for (register int i=0; i<nSpecie; i++)
    {
        W_[i] = this->specieThermo()[i].W();
        href_[i] = 
        (
            this->specieThermo()[i].hc()
          - this->specieThermo()[i].cp
            (
                specie::Pstd,
                specie::Tstd
            )*specie::Tstd
        );
    }

    if (Switch(coeffsDict_.lookup("saveISATtree")))
    {
        this->writeOpt() = IOobject::AUTO_WRITE;
    }


    { // following is needed if old isat table should be read in
      // read in happens at first call to cirxn. has to be in the 
      // right working directory
        if (Pstream::parRun()) //parallel case
        {
            chDir(this->time().path());
        }

        int ci_info[20] = {};
        double ci_rinfo[20] = {};
        int info[100] = {};
        double rinfo[50] = {};
        
        if(Switch(coeffsDict_.lookupOrDefault("constantPressure", false)))
        {
            ci_info[0]=1;
        }
        if(!Switch(this->time().controlDict().lookup("adjustTimeStep")))
        {
            ci_info[1]=1;
        }
        if(Switch(coeffsDict_.lookupOrDefault("externalCKWYP", false)))
        { // if you have put your own ckwyp_ext file in ISAT-CK7 source
            Info << "Using external CKWYP" << endl;
            ci_info[19]=1; 
        }
        {
            ifstream f("isat_1.tab");
            if(f.good()) info[9] = 1; // load an existing table, since it exists. No check if it is valid or not
            f.close();
        }
        {
            ifstream f("streams.in");
            if(!f.good()) // create a streams file here!
            {
                f.close();
                const specie& spec = this->specieThermo()[0];
                const volScalarField& p = mesh.lookupObject<volScalarField>("p");
                std::ofstream streamsFile("streams.in");
                streamsFile << "MODECI         ISAT_DI" << std::endl;
                streamsFile << "STREAM BEGIN DUMMY [MOLE]" << std::endl;
                streamsFile << "    P       " << p[0]/1.01325E+05 << std::endl;
                streamsFile << "    T       300" << std::endl;
                streamsFile << "    " << spec.name() << "      1" << std::endl;
                streamsFile << "STREAM END DUMMY" << std::endl;
                streamsFile.close();
            }
            f.close();
        }
        
        ciparam_( ci_info, ci_rinfo, info, rinfo );
        
        ciinit_( &ncv_, &nfull, &nstrms );
        Info << "ISAT-CK7 initialized" << endl; 

        double deltaT = this->time().deltaT().value();
        double Z0[ncv_];
        double Z1[ncv_];
        int strm=1;
        double dpt[3];

        cistrm_( &strm, &ncv_, Z0, dpt );

        cirxn_(&deltaT, &ncv_, Z0, Z1, dpt); // call once to complete init

        if (Pstream::parRun()) //parallel case
        {
            chDir("..");
        }
    }
    if (ncv_ != nSpecie+1)
    {
        FatalErrorIn("ChemistryModel::ISAT::ISAT(const mesh&)")
            << "The number of species in OpenFOAM does not match "
            << "the number of species in ISAT library:" << nl
            << "Check chem.bin file for errors" << nl
            << abort(FatalError);
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

    int nSpecie = this->nSpecie();
    double Z0[ncv_];  // initial specific mole number [mole/g] array + hs
    double Z1[ncv_];  // specific mole number after deltaT array + hs

    scalar rcTot = 1./sum(c);    

    // Calculate the sensible enthalpy
    typename ChemistryModel::thermoType mixture
    (
        (c[0]*rcTot)*this->specieThermo_[0]
    );

    double densitySI = c[0]*W_[0];

    for (register int i=1; i<nSpecie; i++)
    {
        densitySI += c[i]*W_[i];
        mixture += (c[i]*rcTot)*this->specieThermo_[i];
    }

    scalar hc = 0.;
    for (register int i=0; i<nSpecie; i++)
    {
        Z0[i] = c[i]/densitySI; // kmol/m3 to mole/g
        hc += Z0[i]*href_[i];
    }

    Z0[nSpecie] = (mixture.Ha(p, T)-hc)*1e4; // hs: J/kg to ergs/g

    double dpt[3] = {-1, p*10., -1}; // only p needed as input

    cirxn_(&deltaT, &ncv_, Z0, Z1, dpt);

    for (register int i=0; i<nSpecie; i++)
    {
        c[i] = Z1[i]*densitySI;
    }
    subDeltaT = deltaT;
}

template<class ChemistryModel>
bool Foam::ISAT<ChemistryModel>::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp 
) const
{
    if (Pstream::parRun()) //parallel case
    {
        chDir(this->time().path());
    }
    int nRec = 0;
    cisave_(&nRec);
    if (Pstream::parRun()) //parallel case
    {
        chDir("..");
    }
    return true;
}
// ************************************************************************* //
