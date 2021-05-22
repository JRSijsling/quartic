/***
 *  Random verification of Van Rijnswou's result.
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

SetMemoryLimit(16*2^30);
AttachSpec("../magma/spec");
import "../magma/reconstruction/TernaryForms.m": TernaryToBinary;
load "Ingredients.m";

print "We perform some random sanity checks on Van Rijnswou's results.";

function GL2ToGL3Dual(U);
    // From W to Sym^2 (W^*)
    a,b,c,d := Explode(Eltseq(U));
    return Matrix([[a^2, a*b, b^2], [2*a*c, a*d + b*c, 2*b*d], [c^2, c*d, d^2]]);
end function;

/* Start from a random ternary quartic: */
F := GF(NextPrime(10^20));
//F := Rationals();
R<x1,x2,x3> := PolynomialRing(F, 3);
S<z1,z2> := PolynomialRing(F, 2);
B := 2^3;

for i:=1 to B do
    print "Run", i;
    print "Random quartic form used:";
    f := RandomQuartic(F);
    print f;

    /* Random transformations in degree 2: */
    U := RandomSL2(F);
    print "Random matrix used:";
    print U;

    /* Checking that the homomorphism above fixes the correct form: */
    print "Checking that v_2^2 - v_1 * v_3 is stabilized:";
    print TransformForm(x2^2 - x1*x3, GL2ToGL3Dual(U) : contra := true) eq x2^2 - x1*x3;

    /* Check equivariance of Van Rijnswou's transformation: */
    print "Checking equivariance of transformation:";
    fU := TransformForm(f, GL2ToGL3Dual(U));
    print [ S!g : g in TernaryToBinary(fU) ] eq [ S!TransformForm(g, U) : g in TernaryToBinary(f) ];

    print "";
end for;

exit;
