/***
 *  Verification of Dixmier-Ohno interpolations.
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
import "../magma/invariants/DixmierOhnoInvariants.m": Rho;
import "../magma/reconstruction/TernaryForms.m": TernaryToBinary, QuadricNormalizer, Dehomogenization;
import "../magma/reconstruction/JointCovariantsData.m": S8S4Cov;
import "../magma/reconstruction/JointCovariants.m": IthJointInvariant, JointCovariant;
import "../magma/reconstruction/Interpolations.m": JointInvariantFromDixmierOhno;
load "Ingredients.m";

print "We verify our invariant interpolations.";

/* Parameters for random quartic generation and first sample size: */
coeff_bound := 7;
eps := 5;
char := NextPrime(10^3);

/* Creation of ambient polynomial ring: */
QQ := Rationals();
FF := FiniteField(char);

W := [ 3, 6, 9, 9, 12, 12, 15, 15, 18, 18, 21, 21, 27 ];
R<I03, I06, I09, J09, I12, J12, I15, J15, I18, J18, I21, J21, I27> := PolynomialRing(FF, W);

/* Creation of Hilbert series: */
P<z> := PowerSeriesRing(QQ);
den := (1 - z^3)*(1 - z^6)*(1 - z^9)*(1 - z^12)*(1 - z^15)*(1 - z^18)*(1 - z^27);
num := 1 + z^9 + z^12 + z^15 + 2*z^18 + 3*z^21 + 2*z^24 + 3*z^27 + 4*z^30 +
3*z^33 + 4*z^36 + 4*z^39 + 3*z^42 + 4*z^45 + 3*z^48 + 2*z^51 + 3*z^54 + 2*z^57
+ z^60 + z^63 + z^66 + z^75;
HS := num / den;

/* Full set of invariants used */
names := [ "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "j3", "j4c", "j5e", "j5f", "j6h", "j6i", "j3a", "j4a", "j4b", "j5c", "j5d", "j6e", "j6f", "j6g" ];
degrees := [ 2, 3, 4, 5, 6, 7, 8, 9, 10, 3, 4, 5, 5, 6, 6, 3, 4, 4, 5, 5, 6, 6, 6 ];
indices := [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 74, 79, 85, 86, 95, 96, 75, 77, 78, 83, 84, 92, 93, 94  ];

for counter:=1 to #names do
    name := names[counter];
    degree := degrees[counter];
    index := indices[counter];
    rank := Integers()!Coefficient(HS, 9*degree);
    N := rank + eps;

    print "Verifying", name , "...";
    print "Rank to be attained:", rank;
    n := 0;
    test := true;
    mons := MonomialsOfWeightedDegree(R, 9*degree);
    ev_monss := [];
    print "Generating interpolation points...";
    repeat
        n +:= 1;
        repeat
            f := RandomQuartic(QQ : B := coeff_bound);
            P := ProjectiveSpace(Parent(f));
            rho := Rho(f);
            I := DixmierOhnoInvariants(f);
            I09 := I[3];
            T0 := Matrix([[1/2, 0, -1/2], [0, 1, 0], [1/2, 0, 1/2]]);
            Cs := QuadricCoefficients(TransformForm(rho, T0 : contra := true));
            a := Cs[2,2]; b := Cs[1,2]; d := Cs[1,1];
            delta := (a*d - 1/4*b^2);
        until (not IsSingular(Scheme(P, rho))) and (I09*a*delta ne 0);
        T := QuadricNormalizer(rho);
        f_norm := ChangeRing(f, BaseRing(T));
        f_norm := TransformForm(f_norm, T);
        I_norm := DixmierOhnoInvariants(f_norm);
        b8_norm, b4_norm, b0_norm := Explode(TernaryToBinary(f_norm));
        b0_norm := BaseRing(Parent(b0_norm)) ! b0_norm;

        ev_pol := JointInvariantFromDixmierOhno(name, I);
        if index eq 0 then
            ev_JI := IthJointInvariant(S8S4Cov, [ Dehomogenization(b) : b in [ b4_norm, b8_norm ] ], counter);
            test and:= ( ev_pol / I09^degree ) eq (ev_JI / (b0_norm^degree));
            if n eq 1 then
                print "(Intermediate invariance check:";
                U := RandomSL2(QQ);
                ev_JI_transf := IthJointInvariant(S8S4Cov, [ Dehomogenization(TransformForm(b, U)) : b in [ b4_norm, b8_norm ] ], counter);
                print ev_JI eq ev_JI_transf, ")";
            end if;
        else
            ev_JI := JointCovariant(S8S4Cov, [ Dehomogenization(b) : b in [ b4_norm, b8_norm ] ], index)[1];
            if n eq 1 then
                print "(Intermediate invariance check:";
                U := RandomSL2(QQ);
                ev_JI_transf := JointCovariant(S8S4Cov, [ Dehomogenization(TransformForm(b, U)) : b in [ b4_norm, b8_norm ] ], index)[1];
                print ev_JI eq ev_JI_transf, ")";
            end if;
        end if;

        /* Affine version (not used because of blowup when working over the rationals): */
        //rho_norm := Rho(f_norm);
        //u_norm := Coefficients(Rho)[1];
        //ev_JI := IthJointInvariant(S8S4Cov, [ Dehomogenization(b4_norm), Dehomogenization(b8_norm) ], counter);
        //ev_pol := JointInvariantFromDixmierOhno(name, I_norm);
        //test and:= (ev_pol eq ((5 * u^2)^degree * ev_JI));

        if Min([ Valuation(QQ ! i, char) : i in I ]) ge 0 then
            IFF := [ FF ! i : i in I ];
            Append(~ev_monss, Matrix([ [ FF ! x : x in Eltseq(Evaluate(mon, IFF)) ] : mon in mons ]));
        end if;

        done := false;
        if n ge N then
            M := HorizontalJoin(ev_monss);
            done := (Rank(M) eq rank);
        end if;

    until done;

    print "Sufficiently large set of interpolation points generated";
    print "Test that given expression works:", test;
end for;

exit;
