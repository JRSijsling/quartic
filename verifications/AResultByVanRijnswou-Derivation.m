/***
 *  Rederivation of Van Rijnswou's result.
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

print "We rederive Van Rijnswou's results.";
print "Covariant case:";

// GL2ToGL3:

F := Rationals();
S<z1,z2> := PolynomialRing(F, 2);
R<x1,x2,x3> := PolynomialRing(F, 3);

eS := hom<S -> S | f :-> z2*Derivative(f, z1)>;
hS := hom<S -> S | f :-> z1*Derivative(f, z1) - z2*Derivative(f, z2) >;
fS := hom<S -> S | f :-> z1*Derivative(f, z2)>;
DSs := [* eS, fS, hS *];
eR := hom<R -> R | f :-> 2*x2*Derivative(f, x1) +   x3*Derivative(f, x2)>;
hR := hom<R -> R | f :-> 2*x1*Derivative(f, x1) - 2*x3*Derivative(f, x3) >;
fR := hom<R -> R | f :->   x1*Derivative(f, x2) + 2*x2*Derivative(f, x3)>;
DRs := [* eR, fR, hR *];

// Monomial of weight 8 automatic:
v8 := x1^4;

MonsRd4 := MonomialsOfDegree(R, 4);
MonsSd8 := MonomialsOfDegree(S, 8);
MonsSd4 := MonomialsOfDegree(S, 4);
MonsSd0 := MonomialsOfDegree(S, 0);

// Cannot do this with MonomialsOfWeight (except by shifting in a weird way):
MonsRw4 := [ x1^3*x3, x1^2*x2^2 ];
M4 := Matrix([ [ MonomialCoefficient(fR(mon4), mon) : mon in MonsRd4 ] : mon4 in MonsRw4 ]);
c4 := Eltseq(Basis(Kernel(M4))[1]);
v4 := &+[ c4[i] * MonsRw4[i] : i in [1..#c4] ];

// Cannot do this with MonomialsOfWeight (except by shifting in a weird way):
MonsRw0 := [ x1^2*x3^2, x1*x2^2*x3, x2^4 ];
M0 := Matrix([ [ MonomialCoefficient(fR(mon0), mon) : mon in MonsRd4 ] : mon0 in MonsRw0 ]);
c0 := Eltseq(Basis(Kernel(M0))[1]);
v0 := &+[ c0[i] * MonsRw0[i] : i in [1..#c0] ];

// The ones for binary are easy (note that we have some freedom though):
w8 := z1^8;
w4 := z1^4;
w0 := z1^0;

// Polish:
//v4 := 4*v4;
//v0 := 16*v0;
//v4 := 2*v4;
//v0 := -8*v0;
print "Heighest weight vectors in V:";
print v8;
print v4;
print v0;

T3 := [ ];
T28 := [ ];
T24 := [ ];
T20 := [ ];
while not v8 eq 0 do
    Append(~T3,  [ MonomialCoefficient(v8, mon) : mon in MonsRd4 ]);
    Append(~T28, [ MonomialCoefficient(w8, mon) : mon in MonsSd8 ]);
    v8 := eR(v8);
    w8 := eS(w8);
end while;
while not v4 eq 0 do
    Append(~T3,  [ MonomialCoefficient(v4, mon) : mon in MonsRd4 ]);
    Append(~T24, [ MonomialCoefficient(w4, mon) : mon in MonsSd4 ]);
    v4 := eR(v4);
    w4 := eS(w4);
end while;
while not v0 eq 0 do
    Append(~T3,  [ MonomialCoefficient(v0, mon) : mon in MonsRd4 ]);
    Append(~T20, [ MonomialCoefficient(w0, mon) : mon in MonsSd0 ]);
    v0 := eR(v0);
    w0 := eS(w0);
end while;

// We have the matrices that represent from the right (the images are the
// rows).
// Therefore transposing gives the left action. However, rows rule, so we are
// already happy.
M3 := Matrix(T3);
M28 := Matrix(T28);
M24 := Matrix(T24);
M20 := Matrix(T20);
EqTransfTry :=  M3^(-1)*DiagonalJoin(DiagonalJoin(M28, M24), M20);

print "Matrix representating the equivariant transformation:";
print EqTransfTry;
print "Determinant:", Determinant(EqTransfTry);
print "";

// GL2ToGL3Dual: change matrix (or take derivations contravariantly, as Van
// Rijnswou does by accident)

print "Covariant case:";

F := Rationals();
S<z1,z2> := PolynomialRing(F, 2);
R<x1,x2,x3> := PolynomialRing(F, 3);

eS := hom<S -> S | f :-> z2*Derivative(f, z1)>;
hS := hom<S -> S | f :-> z1*Derivative(f, z1) - z2*Derivative(f, z2) >;
fS := hom<S -> S | f :-> z1*Derivative(f, z2)>;
DSs := [* eS, fS, hS *];
eR := hom<R -> R | f :->   x2*Derivative(f, x1) + 2*x3*Derivative(f, x2)>;
hR := hom<R -> R | f :-> 2*x1*Derivative(f, x1) - 2*x3*Derivative(f, x3) >;
fR := hom<R -> R | f :-> 2*x1*Derivative(f, x2) +   x2*Derivative(f, x3)>;
DRs := [* eR, fR, hR *];

// Monomial of weight 8 automatic:
v8 := x1^4;

MonsRd4 := MonomialsOfDegree(R, 4);
MonsSd8 := MonomialsOfDegree(S, 8);
MonsSd4 := MonomialsOfDegree(S, 4);
MonsSd0 := MonomialsOfDegree(S, 0);

// Cannot do this with MonomialsOfWeight (except by shifting in a weird way):
MonsRw4 := [ x1^3*x3, x1^2*x2^2 ];
M4 := Matrix([ [ MonomialCoefficient(fR(mon4), mon) : mon in MonsRd4 ] : mon4 in MonsRw4 ]);
c4 := Eltseq(Basis(Kernel(M4))[1]);
v4 := &+[ c4[i] * MonsRw4[i] : i in [1..#c4] ];

// Cannot do this with MonomialsOfWeight (except by shifting in a weird way):
MonsRw0 := [ x1^2*x3^2, x1*x2^2*x3, x2^4 ];
M0 := Matrix([ [ MonomialCoefficient(fR(mon0), mon) : mon in MonsRd4 ] : mon0 in MonsRw0 ]);
c0 := Eltseq(Basis(Kernel(M0))[1]);
v0 := &+[ c0[i] * MonsRw0[i] : i in [1..#c0] ];

// The ones for binary are easy (note that we have some freedom though):
w8 := z1^8;
w4 := z1^4;
w0 := z1^0;

// Polish:
v4 := 4*v4;
v0 := 16*v0;
//v4 := 2*v4;
//v0 := -8*v0;
print "Heighest weight vectors in V:";
print v8;
print v4;
print v0;

T3 := [ ];
T28 := [ ];
T24 := [ ];
T20 := [ ];
while not v8 eq 0 do
    Append(~T3,  [ MonomialCoefficient(v8, mon) : mon in MonsRd4 ]);
    Append(~T28, [ MonomialCoefficient(w8, mon) : mon in MonsSd8 ]);
    v8 := eR(v8);
    w8 := eS(w8);
end while;
while not v4 eq 0 do
    Append(~T3,  [ MonomialCoefficient(v4, mon) : mon in MonsRd4 ]);
    Append(~T24, [ MonomialCoefficient(w4, mon) : mon in MonsSd4 ]);
    v4 := eR(v4);
    w4 := eS(w4);
end while;
while not v0 eq 0 do
    Append(~T3,  [ MonomialCoefficient(v0, mon) : mon in MonsRd4 ]);
    Append(~T20, [ MonomialCoefficient(w0, mon) : mon in MonsSd0 ]);
    v0 := eR(v0);
    w0 := eS(w0);
end while;

// We have the matrices that represent from the right (the images are the
// rows).
// Therefore transposing gives the left action. However, rows rule, so we are
// already happy.
M3 := Matrix(T3);
M28 := Matrix(T28);
M24 := Matrix(T24);
M20 := Matrix(T20);
EqTransfTry :=  M3^(-1)*DiagonalJoin(DiagonalJoin(M28, M24), M20);

print "Matrix representating the equivariant transformation:";
print EqTransfTry;
print "Determinant:", Determinant(EqTransfTry);
print "";

exit;
