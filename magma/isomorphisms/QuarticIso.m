/***
 *  Final wrapping intrinsic
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


import "QuarticIsoFF.m": QuarticIsomorphismsFF;
import "QuarticIsoQQ.m": QuarticIsomorphismsQQ;
import "Sutherland.m": SPQIsIsomorphic;


function Normalize33Row(T)

row := Eltseq(Rows(T)[1]);
i0 := Minimum([ i : i in [1..#row] | row[i] ne 0 ]);
return (1/row[i0])*T;

end function;


function Normalize33Column(T)

col := Eltseq(Rows(Transpose(T))[1]);
i0 := Minimum([ i : i in [1..#col] | col[i] ne 0 ]);
return (1/col[i0])*T;

end function;


intrinsic IsIsomorphicQuartic(f1::RngMPolElt, f2::RngMPolElt : geometric := false) -> .
{Tests for the existence of isomorphisms between the ternary quartics f1 and f2, and returns these if they exist.}

K := BaseRing(Parent(f1));
if Type(K) eq FldFin then
    test, Ts := QuarticIsomorphismsFF(f1, f2 : geometric := geometric);
elif Type(K) eq FldRat then
    test, Ts := QuarticIsomorphismsQQ(f1, f2 : geometric := geometric);
else
    test, Ts := SPQIsIsomorphic(f1, f2 : geometric := geometric);
end if;
Ts := [ Normalize33Row(T) : T in Ts ];
return test, Ts;

end intrinsic;


intrinsic AutomorphismGroupQuartic(f::RngMPolElt : geometric := false) -> .
{Finds the automorphism group of the plane quartic curve X as matrices.}

test, Ts := IsIsomorphicQuartic(f, f : geometric := geometric);
return Ts;

end intrinsic;


intrinsic IsIsomorphicQuartic(X1::CrvPln, X2::CrvPln : geometric := false) -> .
{Tests for the existence of isomorphisms between the plane quartic curves X1 and X2, and returns these if they exist.}

f1 := DefiningPolynomial(X1);
f2 := DefiningPolynomial(X2);
assert Degree(f1) eq 4;
assert Degree(f2) eq 4;

test, Ts := IsIsomorphicQuartic(f2, f1 : geometric := geometric);
Ts := [ Normalize33Column(T) : T in Ts ];
return test, Ts;

end intrinsic;


intrinsic AutomorphismGroupQuartic(X::CrvPln : geometric := false) -> .
{Finds the automorphism group of the plane quartic curve X as matrices.}

test, Ts := IsIsomorphicQuartic(X, X : geometric := geometric);
return Ts;

end intrinsic;
