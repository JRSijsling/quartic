/***
 *  Verification of stability claims.
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
import "../magma/reconstruction/TernaryForms.m": BinaryToTernary;
load "Ingredients.m";

print "We verify our statements on stability.";

SetMemoryLimit(16*2^30);

FF := Rationals();
F<f4,f3,f2,f1,f0,g2,g1,g0,h> := FunctionField(FF,9);
R<x1,x2,x3> := PolynomialRing(F,3);
S<z1,z2> := PolynomialRing(F,2);

f := z1^4*(f4*z1^4 + f3*z1^3*z2 + f2*z1^2*z2^2 + f1*z1*z2^3 + f0*z2^4);
g := z1^2*(g2*z1^2 + g1*z1*z2 + g0*z2^2);
h := h;

Phi := R ! BinaryToTernary([f, g, h]);
X := Curve(ProjectiveSpace(R), Phi);
pt := X![0,0,1];

print "Check for non-ordinary singularity:";
print not IsOrdinarySingularity(pt);

exit;
