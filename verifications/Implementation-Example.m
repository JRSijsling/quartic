/***
 *  Reconstruction example used in the paper.
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

AttachSpec("../magma/spec");
AttachSpec("../../hyperelliptic/magma/spec");
SetVerbose("PlaneQuartic", 1);

print "We reconstruct the curve from the paper.";

F := Rationals();
//F := GF(29);

I0 := [F|
    0,
    0,
    0,
    0,
    -133/2687385600,
    0,
    -209/193491763200,
    0,
    2527/23219011584000,
    361/23219011584000,
    -11191/83588441702400000,
    -6137/3482851737600000,
    -2365633/13776537098649600000
];

print "Initial invariants:";
print I0;
f := TernaryQuarticFromDixmierOhnoInvariants(I0);

if F eq Rationals() then
    print "Minimizing defining equation:";
    f := Parent(f) ! MinimizeReducePlaneQuartic(f);
end if;
print "The reconstructed curve is described by the form", f;

I, W := DixmierOhnoInvariants(f);
print "";
print "Test for equality of actual invariants:", I eq I0;
print "";
print "Test for equality of normalized invariants:", WPSNormalize(W, I) eq WPSNormalize(W, I0);



exit;
