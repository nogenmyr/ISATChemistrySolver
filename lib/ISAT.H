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

Class
    Foam::ISAT

Description
    An solver for chemistry using ISAT-CK7
    General note: ISAT-CK7 uses CGS units, so a conversion is done
    before calling. If OpenFOAM is operated using other than SI units,
    this conversion will fail.
    
    ISATCoeffs
    {
        saveISATtree off; //Each processor saves it's generated table if on
        externalCKWYP off; //If you want to use a ckwyp.f file for rates
        constantPressure off; //ISAT can run more efficient if pressure is const
        ISATABabsErr 1e-5; // Abs error for ISAT-tabulation
        ISATABrelErr 0; // Rel error for ISAT-tabulation
        DDASACabsErr 1e-8; // Abs error for chemistry integrator
        DDASACrelErr 1e-9; // Rel error for chemistry integrator
        ISATABsizeMB 500; // Size in megabytes for generated ISAT-table
    }
    

SourceFiles
    ISAT.C

\*---------------------------------------------------------------------------*/

#ifndef ISAT_H
#define ISAT_H

#include "chemistrySolver.H"
#include "IOobject.H"

extern "C" void ciparam_(int*, double*, int*, double*);
extern "C" void ciinit_(int*, int*, int*);
//extern "C" const int*  ciinit_(int*, int*, int*);
extern "C" void cirxn_(double*, int*, double*, double*, double*);
extern "C" void cisave_(int*);
extern "C" void cistrm_(int*, int*, double*, double*);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class ISAT Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class ISAT
:
    public chemistrySolver<ChemistryModel>
{
    // Private data
        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- Mole weights of species
        scalarField W_;
        
        //- Reference enthalpies for species
        scalarField href_;
        
        //- Number of variables for ISAT, nSpecies+1
        mutable int ncv_;


public:

    //- Runtime type information
    TypeName("ISAT");


    // Constructors

        //- Construct from mesh
        ISAT(const fvMesh& mesh);


    //- Destructor
    virtual ~ISAT();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalar& T,
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;


        //- Implements writeObject here to save ISAT table data
        // This is called if writeOpt==IOobject::AUTO_WRITE
        bool writeObject
        (
            IOstream::streamFormat fmt,
		    IOstream::versionNumber ver,
		    IOstream::compressionType cmp 
    	) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ISAT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
