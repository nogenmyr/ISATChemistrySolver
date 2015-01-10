/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2014 Karl-Johan Nogenmyr
-------------------------------------------------------------------------------
License
    Copyright (C) 2014 Karl-Johan Nogenmyr

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

#include "chemistrySolver.H"
#include "chemistryModel.H"

#include "noChemistrySolver.H"
#include "ISAT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistrySolverType(SS, Comp, Thermo)                             \
                                                                              \
    typedef SS<chemistryModel<Comp, Thermo> > SS##Comp##Thermo;               \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        SS##Comp##Thermo,                                                     \
        (#SS"<" + word(Comp::typeName_())                                     \
      + "," + Thermo::typeName() + ">").c_str(),                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        Comp,                                                                 \
        SS##Comp##Thermo,                                                     \
        fvMesh                                                                \
    );


#define makeChemistrySolverTypes(CompChemModel,Thermo)                        \
                                                                              \
    makeChemistrySolverType                                                   \
    (                                                                         \
        ISAT,                                                                 \
        CompChemModel,                                                        \
        Thermo                                                                \
    );                                                                        \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#endif

#include "thermoPhysicsTypes.H"
#include "psiChemistryModel.H"
#include "rhoChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistrySolverTypes(psiChemistryModel, constGasHThermoPhysics);
    makeChemistrySolverTypes(psiChemistryModel, gasHThermoPhysics);
    makeChemistrySolverTypes
    (
        psiChemistryModel,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypes
    (
        psiChemistryModel,
        incompressibleGasHThermoPhysics)
    ;
    makeChemistrySolverTypes(psiChemistryModel, icoPoly8HThermoPhysics);
    makeChemistrySolverTypes(rhoChemistryModel, constGasHThermoPhysics);
    makeChemistrySolverTypes(rhoChemistryModel, gasHThermoPhysics);
    makeChemistrySolverTypes
    (
        rhoChemistryModel,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypes
    (
        rhoChemistryModel,
        incompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypes(rhoChemistryModel, icoPoly8HThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy
    makeChemistrySolverTypes(psiChemistryModel, constGasEThermoPhysics);
    makeChemistrySolverTypes(psiChemistryModel, gasEThermoPhysics);
    makeChemistrySolverTypes
    (
        psiChemistryModel,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypes
    (
        psiChemistryModel,
        incompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypes(psiChemistryModel, icoPoly8EThermoPhysics);
    makeChemistrySolverTypes(rhoChemistryModel, constGasEThermoPhysics);
    makeChemistrySolverTypes(rhoChemistryModel, gasEThermoPhysics);
    makeChemistrySolverTypes
    (
        rhoChemistryModel,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypes
    (
        rhoChemistryModel,
        incompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypes(rhoChemistryModel, icoPoly8EThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
