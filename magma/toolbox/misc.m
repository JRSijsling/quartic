//freeze;

/***
 *  Lifts over the p-adics, or the rationals
 *
 *  Distributed under the terms of the GNU Lesser General Public License (L-GPL)
 *                  http://www.gnu.org/licenses/
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  Copyright 2011-2012, R. Basson & R. Lercier & C. Ritzenthaler & J. Sijsling
 */

 /* Copy of hyperelliptic/magma/toolbox/misc.m */

function MyBaseRing(Rg)
    if IsField(Rg) then
	return BaseField(Rg);
    end if;
    return BaseRing(Rg);
end function;

function MyIsPerfect(Rg)
    if IsField(Rg) then
	return IsPerfect(Rg);
    end if;
    return false;
end function;

function MyIsPrimeField(Rg)
    if IsField(Rg) then
	return IsPrimeField(Rg);
    end if;
    return false;
end function;

forward ChangeBaseRing;
function ChangeBaseRing(elt, F)

    if elt eq 0 then return F!0; end if;

    if MyIsPrimeField(F) then
	if Characteristic(F) eq 0 then
            if Type(Parent(elt)) eq FldRat then
                return F!(Rationals()!elt);
            end if;
	    return F!(IntegerRing()!elt);
	end if;
	return F!(elt);
    end if;

    try
	return F![ChangeBaseRing(x, MyBaseRing(F)) : x in Eltseq(elt)];
    catch e
	nothingdone := 0;
    end try;

    if Rank(F) eq 1 then
	NC := [ ChangeBaseRing(x, MyBaseRing(F)) : x in Coefficients(Numerator(elt)) ];
	N  := &+[NC[i] * F.1^(i-1) : i in [1..#NC]];

	DC := [ ChangeBaseRing(x, MyBaseRing(F)) : x in Coefficients(Denominator(elt)) ];
	D  := &+[DC[i] * F.1^(i-1) : i in [1..#DC]];

	return N/D;
    end if;

    if IsField(F) then
	C, T := CoefficientsAndMonomials(Numerator(elt)); T  := [Exponents(t)  : t in T];
	NC := [ChangeBaseRing(x, MyBaseRing(F)) : x in C];
	NT := [&*[F.i^t[i] : i in [1..#t]]  : t in T];
	N := &+[NC[i] * NT[i] : i in [1..#NT]];

	C, T := CoefficientsAndMonomials(Denominator(elt)); T  := [Exponents(t)  : t in T];
	NC := [ChangeBaseRing(x, MyBaseRing(F)) : x in C];
	NT := [&*[F.i^t[i] : i in [1..#t]]  : t in T];
	D := &+[NC[i] * NT[i] : i in [1..#NT]];

	return N/D;
    end if;



    C, T := CoefficientsAndMonomials(elt); T  := [Exponents(t)  : t in T];
    NC := [ChangeBaseRing(x, MyBaseRing(F)) : x in C];
    NT := [&*[F.i^t[i] : i in [1..#t]]  : t in T];
    N := &+[NC[i] * NT[i] : i in [1..#NT]];

    return N;
end function;

forward LiftRing;
function LiftRing(F : perfect := MyIsPerfect(F), N := 2)

    if MyIsPrimeField(F) then
	if perfect then
	    return pAdicField(Characteristic(F), N);
	else
	    return Rationals();
	end if;
    end if;

    FLp :=  LiftRing(MyBaseRing(F) : perfect := perfect, N := N);
    try
	return ext< FLp | PolynomialRing(FLp)![ChangeBaseRing(c, FLp) : c in Eltseq(DefiningPolynomial(F))]>;
    catch e
	nothingdone := 0;
    end try;

    if IsField(F) then
	if Rank(F) eq 1 then
	    FL := FunctionField(FLp);
	else
	    FL := FunctionField(FLp, Rank(F));
	end if;
    else
	if Rank(F) eq 1 then
	    FL := PolynomialRing(FLp);
	else
	    FL := PolynomialRing(FLp, Rank(F));
	end if;
    end if;

    return FL;

end function;
