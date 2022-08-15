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

    You should have received a_ copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "RedlichKwongPengRobinsonGas.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::RedlichKwongPengRobinsonGas<Specie>::RedlichKwongPengRobinsonGas
(
    const Specie& sp,
    const scalar& Tc,
    const scalar& Vc,
    const scalar& Zc,
    const scalar& Pc,
    const scalar& omega
)
:
    Specie(sp),
    Tc_(Tc),
    Vc_(Vc),
    Zc_(Zc),
    Pc_(Pc),
    omega_(omega)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::RedlichKwongPengRobinsonGas<Specie>::RedlichKwongPengRobinsonGas
(
    const word& name,
    const RedlichKwongPengRobinsonGas& pg
)
:
    Specie(name, pg),
    Tc_(pg.Tc_),
    Vc_(pg.Vc_),
    Zc_(pg.Zc_),
    Pc_(pg.Pc_),
    omega_(pg.omega_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::RedlichKwongPengRobinsonGas <Specie>>
Foam::RedlichKwongPengRobinsonGas<Specie>::clone() const
{
    return autoPtr<RedlichKwongPengRobinsonGas<Specie>>::New(*this);
}


template<class Specie>
inline Foam::autoPtr<Foam::RedlichKwongPengRobinsonGas<Specie>>
Foam::RedlichKwongPengRobinsonGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<RedlichKwongPengRobinsonGas<Specie>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);
    return p/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::H(scalar p, scalar T) const
{
    const scalar v = this->Z(p, T)*RR*T/p;

    return
        (p*v
        - RR*T
        + 1/(b_*d1md2_)*log((v + d1_*b_)/(v + d2_*b_))
        *(a_*this->alpha(T) - T*this->ap(T) ));
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::Cp(scalar p, scalar T) const
{
    const scalar A = this->A(p, T);
    const scalar B = this->B(p, T);

    const scalar Z = this->Z(p, T);

    const scalar ap = this->ap(T);
    const scalar app = this->app(T);

    const scalar M = (sqr(Z) + d1pd2_*B*Z + d1td2_*sqr(B))/(Z - B);
    const scalar N = ap*B/(b_*RR);

    return
    (
        app*(T/((d1md2_)*b_))*log((Z + d1_*B)/(Z + d2_*B))
      + RR*sqr(M - N)/(sqr(M) - A*(2*Z + d1pd2_*B))
      - RR
    )/this->W();
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::E(scalar p, scalar T) const
{
    const scalar v = this->Z(p, T)*RR*T/p;

    return
        1/(b_*d1md2_)*log((v + d1_*b_)/(v + d2_*b_))
        *(a_*this->alpha(T) - T*this->ap(T) );
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::Cv(scalar p, scalar T) const
{

    const scalar B = this->B(p,T);

    const scalar Z = this->Z(p, T);

    const scalar app = this->app(T);

    return
    (
        app*(T/((d1md2_)*b_))*log((Z + d1_*B)/(Z  + d2_*B))
      - RR
    )/this->W();
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::S
(
    scalar p,
    scalar T
) const
{
    const scalar v = this->Z(p, T)*RR*T/p;

    return
        RR*log(1- b_/v)
        - 1/(b_*d1md2_)*log((v + d1_*b_)/(v + d2_*b_))
        *T*this->ap(T);
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);

    return 1.0/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::Z
(
    scalar p,
    scalar T
) const
{

    const scalar A = this->A(p,T);
    const scalar B = this->B(p,T);

    const scalar a2 = B*(d1pd2_ - 1)- 1;
    const scalar a1 = A + B*(d1td2_*B - d1pd2_*(B + 1));
    const scalar a0 = -B*(d1td2_*(sqr(B) + B) + A);

    const scalar Q = (3*a1 - a2*a2)/9.0;
    const scalar Rl = (9*a2*a1 - 27*a0 - 2*a2*a2*a2)/54.0;

    const scalar Q3 = Q*Q*Q;
    const scalar D = Q3 + Rl*Rl;

    scalar root = -1;

    if (D <= 0)
    {
        const scalar th = ::acos(Rl/sqrt(-Q3));
        const scalar qm = 2*sqrt(-Q);
        const scalar r1 = qm*cos(th/3.0) - a2/3.0;
        const scalar r2 =
            qm*cos((th + 2*constant::mathematical::pi)/3.0) - a2/3.0;
        const scalar r3 =
            qm*cos((th + 4*constant::mathematical::pi)/3.0) - a2/3.0;

        root = max(r1, max(r2, r3));
    }
    else
    {
        // One root is real
        const scalar D05 = sqrt(D);
        const scalar S = cbrt(Rl + D05);
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -cbrt(mag(Rl - D05));
        }
        else
        {
            Tl = cbrt(Rl - D05);
        }

        root = S + Tl - a2/3.0;
    }
    return root;
}


template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    const scalar A = this->A(p,T);
    const scalar B = this->B(p,T);

    const scalar Z = this->Z(p, T);

    const scalar ap = this->ap(T);
    const scalar M = (sqr(Z) + d1pd2_*B*Z + d1td2_*sqr(B))/(Z - B);
    const scalar N = ap*B/(b_*RR);

    return this->R()*sqr(M - N)/(sqr(M) - A*(2*Z + d1pd2_*B));
}

template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::alpha
(
    scalar T
) const
{
    return pow((3.0/(2 + T/Tc_)),kappa_); 
}

template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::ap
(
    scalar T
) const
{
    return -a_*pow(3.0,kappa_)*kappa_/(Tc_*pow((2 + T/Tc_),(kappa_ + 1)));
}

template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::app
(
    scalar T
) const
{
    return a_*(pow(3.0,kappa_)*kappa_*(kappa_ + 1)
                                /(sqr(Tc_)*pow((2 + T/Tc_),(kappa_ + 2))));
}

template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::A
(   
    scalar p, 
    scalar T
) const
{
    return a_*this->alpha(T)*p/sqr(RR*T);
}

template<class Specie>
inline Foam::scalar Foam::RedlichKwongPengRobinsonGas<Specie>::B
(
    scalar p, 
    scalar T
) const
{
    return b_*p/(RR*T);
}
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::RedlichKwongPengRobinsonGas<Specie>::operator+=
(
    const RedlichKwongPengRobinsonGas<Specie>& pg
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pg);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pg.Y()/this->Y();

        Tc_ = Y1*Tc_ + Y2*pg.Tc_;
        Vc_ = Y1*Vc_ + Y2*pg.Vc_;
        Zc_ = Y1*Zc_ + Y2*pg.Zc_;
        Pc_ = RR*Zc_*Tc_/Vc_;
        omega_ = Y1*omega_ + Y2*pg.omega_;
    }
}


template<class Specie>
inline void Foam::RedlichKwongPengRobinsonGas<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::RedlichKwongPengRobinsonGas<Specie> Foam::operator+
(
    const RedlichKwongPengRobinsonGas<Specie>& pg1,
    const RedlichKwongPengRobinsonGas<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return RedlichKwongPengRobinsonGas<Specie>
        (
            sp,
            pg1.Tc_,
            pg1.Vc_,
            pg1.Zc_,
            pg1.Pc_,
            pg1.omega_
        );
    }
    else
    {
        const scalar Y1 = pg1.Y()/sp.Y();
        const scalar Y2 = pg2.Y()/sp.Y();

        const scalar Tc = Y1*pg1.Tc_ + Y2*pg2.Tc_;
        const scalar Vc = Y1*pg1.Vc_ + Y2*pg2.Vc_;
        const scalar Zc = Y1*pg1.Zc_ + Y2*pg2.Zc_;

        return RedlichKwongPengRobinsonGas<Specie>
        (
            sp,
            Tc,
            Vc,
            Zc,
            RR*Zc*Tc/Vc,
            Y1*pg1.omega_ + Y2*pg2.omega_
        );
    }
}


template<class Specie>
Foam::RedlichKwongPengRobinsonGas<Specie> Foam::operator*
(
    const scalar s,
    const RedlichKwongPengRobinsonGas<Specie>& pg
)
{
    return RedlichKwongPengRobinsonGas<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.Tc_,
        pg.Vc_,
        pg.Zc_,
        pg.Pc_,
        pg.omega_
    );
}


template<class Specie>
Foam::RedlichKwongPengRobinsonGas<Specie> Foam::operator==
(
    const RedlichKwongPengRobinsonGas<Specie>& pg1,
    const RedlichKwongPengRobinsonGas<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
     == static_cast<const Specie&>(pg2)
    );

    const scalar Y1 = pg1.Y()/sp.Y();
    const scalar Y2 = pg2.Y()/sp.Y();

    const scalar Tc = Y2*pg2.Tc_ - Y1*pg1.Tc_;
    const scalar Vc = Y2*pg2.Vc_ - Y1*pg1.Vc_;
    const scalar Zc = Y2*pg2.Zc_ - Y1*pg1.Zc_;

    return RedlichKwongPengRobinsonGas<Specie>
    (
        sp,
        Tc,
        Vc,
        Zc,
        RR*Zc*Tc/Vc,
        Y2*pg2.omega_ - Y1*pg1.omega_
    );
}


// ************************************************************************* //
