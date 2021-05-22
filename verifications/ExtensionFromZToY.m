/***
 *  Verification of rationality and integrality properties when using T_int.
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

print "We check integrality properties of a specific normalizing matrix.";

R<a00,a01,a02,a11,a12,a22,a400,a310,a301,a220,a211,a202,a130,a121,a112,a103,a040,a031,a022,a013,a004> := PolynomialRing(Rationals(), 21);
F<a00,a01,a02,a11,a12,a22,a400,a310,a301,a220,a211,a202,a130,a121,a112,a103,a040,a031,a022,a013,a004> := FieldOfFractions(R);
S<x1,x2,x3> := PolynomialRing(F, 3);

as := GeneratorsSequence(F);
mons2 := MonomialsOfDegree(S, 2);
mons4 := MonomialsOfDegree(S, 4);

Phi := &+[ as[i] * mons4[i-6] : i in [7..21] ];
Rho := &+[ as[i] * mons2[i] : i in [1..6] ];

print "Checking polynomiality of square of determinant of normalizer:";
N := IntegralQuadricNormalizer(Rho);
print Determinant(N)^2 in R;
//print R ! (Determinant(N)^2);

print "Checking polynomiality of b0:";
Phi := ChangeRing(Phi, Parent(N[1,1]));
Phi_norm := TransformForm(Phi, N);
b8, b4, b0 := Explode(TernaryToBinary(Phi_norm));
b0 := Coefficients(b0)[1];
print b0 in R;
//print R ! b0;

exit;
