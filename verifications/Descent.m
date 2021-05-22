/***
 *  Random verification of the fundamental lemma.
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

F := Rationals();
R<t> := PolynomialRing(F);
F<r> := NumberField(t^2 - 2);
R<x1,x2,x3> := PolynomialRing(F, 3);

f0 := x1^4 + x1^2*x2*x3 + x1*x2*x3^2 + x1*x3^3 + x2^4 - x3^4;
T := Matrix(F, [[1, r, 1], [2, 0, 1], [r, 1, 0]]);
T := Matrix(F, [[1, r, 0], [1, 1, 0], [0, 0, 1]]);
Ts := Matrix(F, [[1, -r, 0], [1, 1, 0], [0, 0, 1]]);
f := TransformForm(f0, T);
fs := TransformForm(f0, Ts);

print "Verification of calculations in the section on descent:";

print "The ternary quartic f:";
print f;
print "The conjugate of the ternary quartic f:";
print fs;

print "Check that isomorphism works:";
print fs eq TransformForm(f, T^(-1)*Ts);
print "Matrix reprsenting the corresponding cocycle:";
print T^(-1)*Ts;
print "Check that we get back to a ternary quartic over QQ:";
print TransformForm(f, T^(-1));

exit;
