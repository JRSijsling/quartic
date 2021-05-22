/***
 *  Random verification of divisibility properties in main theorem
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
import "../magma/invariants/DixmierOhnoInvariants.m": Rho, ContravariantSigmaAndPsi, DifferentialOperation, JOperation11;
import "../magma/reconstruction/TernaryForms.m": TernaryToBinary;
load "Ingredients.m";

print "We check that the generic divisility properties by I12 of the joint invariants hold.";

F := FiniteField(NextPrime(53));
//F := Rationals();

repeat
    print "Random quartic form used:";
    repeat
        f := RandomQuartic(F);
        P := ProjectiveSpace(Parent(f));
        rho := Rho(f);
    until IsSingular(Scheme(P, rho));
    print f;

    N := IntegralQuadricNormalizer(rho);
    f := ChangeRing(f, Parent(N[1,1]));
    f_norm := TransformForm(f, N);
    b8, b4, b0 := Explode(TernaryToBinary(f_norm));
    I := DixmierOhnoInvariants(f);
    I03, I06, I09, J09, I12, J12, I15, J15, I18, J18, I21, J21, I27 := Explode(I);
    done := false;
    if b0 ne 0 then
        done := true;
        "Test done! I12 is zero while b0 is not.";
    end if;
until done;

exit;
