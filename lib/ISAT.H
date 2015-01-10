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

SourceFiles
    ISAT.C

\*---------------------------------------------------------------------------*/

#ifndef ISAT_H
#define ISAT_H

#include "chemistrySolver.H"
extern "C" void ciinit_(int*, int*, int*);
extern "C" void cirxn_(double*, int*, double*, double*, double*);
extern "C" void ciconv_
(
    int*, int*, int*, int*, int*, double*, double*, double*, 
    double*, int*, int*, int*, int*, double*, double*, double*, double*
);

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
        
        scalarField W_;

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
