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
import "../magma/invariants/DixmierOhnoInvariants.m": Rho, Tau, DifferentialOperation, JOperation11;
import "../magma/reconstruction/TernaryForms.m": QuadricNormalizer, TernaryToBinary;
load "Ingredients.m";

print "We give a random check of the equalities in the fundamental lemma.";

F := FiniteField(53);
F := Rationals();
B := 2^3;

for i:=1 to B do
    print "Run", i;
    print "Random quartic form used:";
    repeat
        f := RandomQuartic(F);
        P := ProjectiveSpace(Parent(f));
        rho := Rho(f);
    until not IsSingular(Scheme(P, rho));
    print f;

    T := QuadricNormalizer(rho);
    detT := Determinant(T);
    R<x1, x2, x3> := PolynomialRing(BaseRing(T), 3);
    f := R ! f;
    rho := R ! rho;
    f_tilde := TransformForm(f, T);
    rhoT := TransformForm(rho, T : contra := true);
    print "Rho (a, v T):";
    print rhoT;

    print "Verifying equality for rho tilde:";
    u := MonomialCoefficient(rhoT, x2^2);
    rho_tilde := Rho(f_tilde);
    print rho_tilde eq (detT^6 * u * (x2^2 - x1*x3));
    print "Verifying equality for tau tilde:";
    tau_tilde := Tau(f_tilde);
    print tau_tilde eq (detT^6 * DifferentialOperation(rhoT, f_tilde) /12);
    print "Verifying first equality for I09 tilde:";
    I_tilde := DixmierOhnoInvariants(f_tilde);
    I09_tilde := I_tilde[3];
    print I09_tilde eq (detT^6 * JOperation11(tau_tilde, rhoT));
    print "Verifying second equality for I09 tilde:";
    a202_tilde := MonomialCoefficient(f_tilde, x1^2*x3^2);
    a121_tilde := MonomialCoefficient(f_tilde, x1*x2^2*x3);
    a040_tilde := MonomialCoefficient(f_tilde, x2^4);
    print I09_tilde eq (detT^12 * u^2 * (a202_tilde/6 - a121_tilde/6 + a040_tilde));
    print "Verifying first equality for I09:";
    I := DixmierOhnoInvariants(f);
    I09 := I[3];
    print I09 eq (u^2 * (a202_tilde/6 - a121_tilde/6 + a040_tilde));
    print "Verifying equality for b0 l^* T:";
    b0lT := TernaryToBinary(f_tilde)[3];
    print b0lT eq (a202_tilde/30 - a121_tilde/30 + a040_tilde/5);
    print "Verifying second equality for I09:";
    print I09 eq 5 * u^2 * b0lT;
    print "Verifying equality between u and I_12:";
    I12 := I[5];
    print -(u^3 / 4) eq (I12 / detT^2);
    u := MonomialCoefficient(rho_tilde, x2^2);
    print "Lemma: verifying whether I09 equals 5 u^2 b0l on Z:";
    print I09_tilde eq (5 * u^2 * b0lT);
    print "Lemma: verifying whether (I09 / b0lT)^3 equals 2000 (I12^2 / det(T)^2)^2 on Y:";
    print (I09 / b0lT)^3 eq (2000 * (I12 / detT^2)^2);
    print "";
end for;

exit;
