/***
 *  Various algorithms for generating random forms, plus integral normalization
 *  matrix
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

/* Generating random objects: */

function RandomDomain(F : B := 5)
    if IsFinite(F) then
        D := F;
    else
        D := [ -B..B ];
    end if;
    return D;
end function;

function RandomNonZero(D : B := 5)
    repeat
        r := Random(D);
    until r ne 0;
    return r;
end function;

function RandomConic(F : B := 5)
    R<x1, x2, x3> := PolynomialRing(F, 3);
    D := RandomDomain(F : B := B);
    return &+[ Random(D) * mon : mon in MonomialsOfDegree(R, 2) ];
end function;

function RandomQuartic(F : B := 5)
    R<x1, x2, x3> := PolynomialRing(F, 3);
    D := RandomDomain(F : B := B);
    return &+[ Random(D) * mon : mon in MonomialsOfDegree(R, 4) ];
end function;

function RandomGL2(F : B := 5)
    D := RandomDomain(F : B := B);
    repeat
        T := Matrix(F, 2, 2, [ Random(D) : i in [1..2*2] ]);
    until Determinant(T) ne 0;
    return T;
end function;

function RandomSL2(F : B := 5)
    D := RandomDomain(F : B := B);
    Diag := [ F!RandomNonZero(D) : i in [1..1] ];
    Append(~Diag, (&*Diag)^(-1));
    T := RandomGL2(F);
    return T * DiagonalMatrix(Diag) * T^(-1);
end function;

function RandomGL3(F: B := 5)
    D := RandomDomain(F : B := B);
    repeat
        T := Matrix(F, 3, 3, [ Random(D) : i in [1..3*3] ]);
    until Determinant(T) ne 0;
    return T;
end function;

function RandomSL3(F : B := B)
    D := RandomDomain(F : B := B);
    Diag := [ F!RandomNonZero(D) : i in [1..2] ];
    Append(~Diag, (&*Diag)^(-1));
    T := RandomGL3(F);
    return T * DiagonalMatrix(Diag) * T^(-1);
end function;

/* Unsophisticated version of an integral normalizing matrix: */
/* (A better one is in my notes, but we will not need it.) */

function QuadricCoefficients(Q)
    S := Parent(Q);
    n := #Names(S);
    return Matrix([ [ MonomialCoefficient(Q, S.i*S.j) : j in [1..n] ] : i in [1..n] ]);
end function;

function IntegralQuadricNormalizer(Q)
    R := BaseRing(Parent(Q));
    T0 := Matrix([[1/2, 0, -1/2], [0, 1, 0], [1/2, 0, 1/2]]);
    Cs := QuadricCoefficients(TransformForm(Q, T0 : contra := true));
    a := Cs[2,2]; b := Cs[1,2]; c := Cs[2,3]; d := Cs[1,1]; e := Cs[1,3]; f := Cs[3,3];
    delta := (a*d - 1/4*b^2);
    N := (a*e - 1/2*b*c);
    Delta := (a*d*f - 1/4*a*e^2 - 1/4*b^2*f + 1/4*b*c*e - 1/4*c^2*d);
    test, r := IsSquare(-a*Delta);
    if not test then
        P<u> := PolynomialRing(R);
	r := FieldOfFractions(quo<P | u^2 + a*Delta>).1;
    end if;
    // Weight 0 version that restricts to identity:
    T1 := a*Matrix([[1, b/(2*a), 0], [0, 1, 0], [0, c/(2*a), 1]]);
    T2 := delta*Matrix([[1, 0, 0], [0, 1, 0], [N/(2*delta), 0, 1]]);
    // Inbetween these two steps we could use a diagonal matrix instead of the
    // multiplications above: this would be more efficient.
    T3 := delta*Matrix([[1, 0, 1], [0, 1, 0], [r/delta, 0, -r/delta]]);
    T4 := a^2*Matrix([[-delta/a^2, 0, 0], [0, 1, 0], [0, 0, 1]]);
    return T0*T1*T2*T3*T4;
end function;
