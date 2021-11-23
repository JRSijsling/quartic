/***
 *  Generating random examples of reconstruction
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

SetVerbose("QuarticRec", 1);
SetVerbose("PlaneQuartic", 1);
SetVerbose("ClusterReduction", 1);

/* Start from a random ternary quartic */
F := GF(23);/* Rationals(); */  B := 2^10; Domain := [-B..B];
//R<x> := PolynomialRing(Rationals()); F<r> := NumberField(x^2 - 2);  B := 2^3; Domain := [-B..B];
//F := GF(NextPrime(10^5));  Domain := F;
R<x1, x2, x3> := PolynomialRing(F, 3);
repeat
    f := &+[ (F ! Random(Domain))*mon : mon in MonomialsOfDegree(R, 4) ];
    DOf, DOWght := DixmierOhnoInvariants(f);
until DOf[#DOf] ne 0;

/* Display curve and normalize invariants: */
print "";
print "Start from f =", f;
print "";
print "Its invariants are", DOf;
DOf_norm := WPSNormalize(DOWght, DOf);
/* TODO: Introducing a single incorrect element still gives output */
//DOf_norm[2] := 10;
//print DOf_norm;

/* Construct another quartic with equivalent invariants */
g, aut, twists := TernaryQuarticFromDixmierOhnoInvariants(DOf);
_<x1,x2,x3> := Parent(g);
print "";
print "The reconstructed curve is g =", g;
print "Automorphism group", aut;

/* Compute its invariants and normalize */
DOg, _ := DixmierOhnoInvariants(g);
ChangeUniverse(~DOg, F);
DOg_norm := WPSNormalize(DOWght, DOg);

/* Check everything's fine */
print "";
print "Test for equality of invariants:", DOf eq DOg;
print "";
print "Test for equality of normalized invariants:", DOf_norm eq DOg_norm;
print DOg_norm;
