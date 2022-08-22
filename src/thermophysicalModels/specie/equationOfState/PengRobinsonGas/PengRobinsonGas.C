/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
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

#include "PengRobinsonGas.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::PengRobinsonGas<Specie>::PengRobinsonGas
(
    const dictionary& dict
)
:
    Specie(dict),
    Tc_(dict.subDict("equationOfState").lookup<scalar>("Tc")),
    Vc_(dict.subDict("equationOfState").lookup<scalar>("Vc")),
    Zc_(1.0),
    Pc_(dict.subDict("equationOfState").lookup<scalar>("Pc")),
    omega_(dict.subDict("equationOfState").lookup<scalar>("omega"))
{
    Zc_ = Pc_*Vc_/(RR*Tc_);
    
    // EoS coefficients
    a_coeff_ = 0.45724;
    b_coeff_ = 0.07780;
    a_ = a_coeff_*sqr(RR*Tc_)/Pc_;
    b_ = b_coeff_*RR*Tc_/Pc_;
    
    d1_ = 1 + sqrt(2.0);
    d2_ = 1 - sqrt(2.0);

    d1pd2_ = d1_ + d2_; // == u parameter in EoS
    d1td2_ = d1_*d2_; // == w parameter in EoS
    d1md2_ = d1_ - d2_;

    kappa_ = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::PengRobinsonGas<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PengRobinsonGas<Specie>& pg
)
{
    pg.write(os);
    return os;
}


// ************************************************************************* //
