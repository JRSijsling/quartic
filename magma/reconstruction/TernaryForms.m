/***
 *  Various algorithms for ternary forms
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
 *  Copyright 2016 R. Lercier, C. Ritzenthaler & J.R. Sijsling
 */

 /***
 * Exported intrinsics.
 *
 * intrinsic TransformForm(f::RngMPolElt, T::AlgMatElt :
 *     co := true, contra := false) -> .
 *
 ********************************************************************/

/* Current transformation algorithm: */

intrinsic TransformForm(f::RngMPolElt, T::AlgMatElt :
    co := true, contra := false) -> .
{Transforms the form f by the matrix T.}
    R := Parent(f);
    vars := Matrix([ [ mon ] : mon in MonomialsOfDegree(R, 1) ]);
    if (not co) or contra then
        return Evaluate(f, Eltseq(ChangeRing(Transpose(T)^(-1), R) * vars));
    end if;
    return Evaluate(f, Eltseq(ChangeRing(T, R) * vars));
end intrinsic;

/* From form to polynomials and back: */

function Dehomogenization(b);
    S := PolynomialRing(BaseRing(Parent(b))); t := S.1;
    return S ! Evaluate(b, [t, 1]);
end function;

function Homogenization(f : degree := Degree(f));
    R := PolynomialRing(BaseRing(Parent(f)), 2); z1 := R.1; z2 := R.2;
    return R ! (z2^degree * Evaluate(f, z1/z2));
end function;

/* Normalizing quadrics: */

function QuadricCoefficients(Q)
    S := Parent(Q);
    n := #Names(S);
    return Matrix([ [ MonomialCoefficient(Q, S.i*S.j) : j in [1..n] ] : i in [1..n] ]);
end function;

function QuadricNormalizer(Q)
    // Take care when using this function: it is to be applied to the form
    // itself to normalize the covariant.
    R := BaseRing(Parent(Q));
    // This first step is to ensure that it works on forms already normalized:
    T0 := Matrix([[1/2, 0, -1/2], [0, 1, 0], [1/2, 0, 1/2]]);
    Cs := QuadricCoefficients(TransformForm(Q, T0 : contra := true));
    a := Cs[2,2]; b := Cs[1,2]; c := Cs[2,3]; d := Cs[1,1]; e := Cs[1,3]; f := Cs[3,3];
    delta := (a*d - 1/4*b^2);
    // Do we have any trouble ?
    if not ( IsUnit(a) and IsUnit(delta) ) then
        return Matrix([ [0, 0, 0], [0, 0, 0],  [0, 0, 0] ]);
    end if;
    N := (a*e - 1/2*b*c);
    Delta := (a*d*f - 1/4*a*e^2 - 1/4*b^2*f + 1/4*b*c*e - 1/4*c^2*d);
    test, r := IsSquare(-a*Delta);
    if not test then
        P := PolynomialRing(R); u := P.1;
        r := FieldOfFractions(quo<P | u^2 + a*Delta>).1;
    end if;
    T1 := Matrix([[1, b/(2*a), 0], [0, 1, 0], [0, c/(2*a), 1]]);
    T2 := Matrix([[1, 0, 0], [0, 1, 0], [N/(2*delta), 0, 1]]);
    T3 := Matrix([[1, 0, 1], [0, 1, 0], [r/delta, 0, -r/delta]]);
    T4 := Matrix([[-delta/a^2, 0, 0], [0, 1, 0], [0, 0, 1]]);
    return T0*T1*T2*T3*T4;
end function;

/* From ternary to binary and back */

function TernaryToBinary(f);
    R := Parent(f);
    F := BaseRing(R);
    S := PolynomialRing(F, 2);
    T := Matrix(F, [
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, -3/14, 0, 0, 0, 0, 0],
        [0, 0, 4, 0, 0, 0, 0, 0, 0, 1/7, 0, 0, 0, 0, 0],
        [0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -2/7, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2/7, 0, 0, 1/30],
        [0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 6/7, 0, 0, 0, 0],
        [0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, -1/7, 0, 0, -1/30],
        [0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -2/7, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -3/14, 0],
        [0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 24/7, 0, 0, 1/5],
        [0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 6/7, 0, 0],
        [0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 1/7, 0],
        [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
    ]);
    monsR := MonomialsOfDegree(R, 4);
    v := Matrix(F, [ [ MonomialCoefficient(f, mon) : mon in monsR ] ]);
    vT := v*T;
    monsS8 := MonomialsOfDegree(S, 8);
    monsS4 := MonomialsOfDegree(S, 4);
    monsS0 := MonomialsOfDegree(S, 0);
    b8 := &+[ vT[1,i]      * monsS8[i] : i in [1..9] ];
    b4 := &+[ vT[1,i + 9]  * monsS4[i] : i in [1..5] ];
    b0 := &+[ vT[1,i + 14] * monsS0[i] : i in [1..1] ];
    return [b8, b4, b0];
end function;

function BinaryToTernary(bs);
    b8, b4, b0 := Explode(bs);
    S := Parent(b8);
    F := BaseRing(S);
    R := PolynomialRing(F, 3);
    v8 := [ MonomialCoefficient(b8, monb) : monb in MonomialsOfDegree(S, 8) ];
    v4 := [ MonomialCoefficient(b4, monb) : monb in MonomialsOfDegree(S, 4) ];
    v0 := [ MonomialCoefficient(b0, monb) : monb in MonomialsOfDegree(S, 0) ];
    v := Matrix(F, [ v8 cat v4 cat v0 ]);
    Tinv := Matrix(F, [
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, -3/14, 0, 0, 0, 0, 0],
        [0, 0, 4, 0, 0, 0, 0, 0, 0, 1/7, 0, 0, 0, 0, 0],
        [0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -2/7, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2/7, 0, 0, 1/30],
        [0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 6/7, 0, 0, 0, 0],
        [0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, -1/7, 0, 0, -1/30],
        [0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -2/7, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -3/14, 0],
        [0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 24/7, 0, 0, 1/5],
        [0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 6/7, 0, 0],
        [0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 1/7, 0],
        [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
    ])^(-1);
    w := v * Tinv;
    monst := MonomialsOfDegree(R, 4);
    f := &+[ w[1, i] * monst[i] : i in [1..15] ];
    return f;
end function;

/* Conjugating objects by Galois automorphisms: */

function ConjugateForm(sigma, f);
    return &+[ sigma(MonomialCoefficient(f, mon)) * mon : mon in Monomials(f) ];
end function;

function ConjugateMatrix(sigma, M);
    return Matrix([ [ sigma(c) : c in Eltseq(row) ] : row in Rows(M) ]);
end function;
