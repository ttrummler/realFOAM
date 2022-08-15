/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "psiThermo.H"
#include "makeThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "PengRobinsonGas.H"
#include "SoaveRedlichKwongGas.H"
#include "RedlichKwongPengRobinsonGas.H"
#include "hConstThermo.H"
#include "eConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

#include "hePsiThermo.H"
#include "pureMixture.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * * Enthalpy-based * * * * * * * * * * * * * */

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);
// PR
makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    PengRobinsonGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    PengRobinsonGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    janafThermo,
    PengRobinsonGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    PengRobinsonGas,
    specie
);

// SRK
makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    SoaveRedlichKwongGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    SoaveRedlichKwongGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    janafThermo,
    SoaveRedlichKwongGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    SoaveRedlichKwongGas,
    specie
);

// RKPR

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    RedlichKwongPengRobinsonGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    RedlichKwongPengRobinsonGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    janafThermo,
    RedlichKwongPengRobinsonGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    RedlichKwongPengRobinsonGas,
    specie
);
/* * * * * * * * * * * * * * Internal-energy-based * * * * * * * * * * * * * */

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    eConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    eConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    perfectGas,
    specie
);

// PR
makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    PengRobinsonGas,
    specie
);

// SRK
makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    SoaveRedlichKwongGas,
    specie
);

// RKPR
makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    RedlichKwongPengRobinsonGas,
    specie
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
