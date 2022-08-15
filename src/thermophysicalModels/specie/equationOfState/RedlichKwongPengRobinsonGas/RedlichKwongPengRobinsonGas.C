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

#include "RedlichKwongPengRobinsonGas.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::RedlichKwongPengRobinsonGas<Specie>::RedlichKwongPengRobinsonGas
(
    const dictionary& dict
)
:
    Specie(dict),
    Tc_(dict.subDict("equationOfState").get<scalar>("Tc")),
    Vc_(dict.subDict("equationOfState").get<scalar>("Vc")),
    Zc_(1.0),
    Pc_(dict.subDict("equationOfState").get<scalar>("Pc")),
    omega_(dict.subDict("equationOfState").get<scalar>("omega"))
{
    Zc_ = Pc_*Vc_/(RR*Tc_);
    
    // For calculation of EoS coefficients (see Kim et al.,2012)
    scalar Zc_here = Zc_*1e-3; 
    scalar c_z = 1.168;
    scalar d1 = 0.428363;
    scalar d2 = 18.496215;
    scalar d3 = 0.338426;
    scalar d4 = 0.660000;
    scalar d5 = 789.723105;
    scalar d6 = 2.512392;

    d1_ = (d1 + d2*exp(d4*log(d3 - c_z*Zc_here))
                  + d5*exp(d6*log(d3 - c_z*Zc_here)));
    
    d2_ = (1 - d1_)/(1 + d1_);

    d1pd2_ = d1_ + d2_; // == u
    d1td2_ = d1_*d2_; // == w
    d1md2_ = d1_ - d2_;

    scalar d_rkpr = (1 + sqr(d1_))/(1 + d1_);
    scalar y_rkpr = 1 + pow((2*(1 + d1_)),(1./3.)) 
                    + pow((4/(1 + d1_)),(1./3.));

    a_coeff_ = ((3*sqr(y_rkpr) + 3*y_rkpr*d_rkpr 
           + sqr(d_rkpr) + d_rkpr - 1)
                   /sqr(3*y_rkpr + d_rkpr - 1));

    b_coeff_ = 1/(3*y_rkpr + d_rkpr - 1);

    a_ = a_coeff_*sqr(RR*Tc_)/Pc_;
    b_ = b_coeff_*RR*Tc_/Pc_;

    scalar a1 =  66.125;
    scalar a0 = -23.359;
    scalar b1 = -40.594;
    scalar b0 =  16.855;
    scalar c1 =  5.27345;
    scalar c0 = -0.25826;
    
    scalar coeff1 = c_z*Zc_here*a1 + a0;
    scalar coeff2 = c_z*Zc_here*b1 + b0;
    scalar coeff3 = c_z*Zc_here*c1 + c0;
    
    // kappa (or c_alpha)
    // NOTE that the formumlation differs from PR and SRK
    kappa_ = coeff1*sqr(omega_) + coeff2*omega_ + coeff3;

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::RedlichKwongPengRobinsonGas<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const RedlichKwongPengRobinsonGas<Specie>& pg
)
{
    pg.write(os);
    return os;
}


// ************************************************************************* //
