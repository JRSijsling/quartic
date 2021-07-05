/***
 *  Computes the Dixmier-Ohno invariants of a ternary quartic in char. 2
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
 *  Copyright:
 *  2020 R. Lercier, C. Ritzenthaler & J.R. Sijsling
 */

 import "DixmierOhnoInvariants.m" : DerivativeSequence,  PowerDerivative, DifferentialOperation, JOperation11, JOperation22, JOperation30, JOperation03, CovariantHessian,    ContravariantSigmaAndPsi, DixmierInvariant;

import "../toolbox/misc.m" : ChangeBaseRing, LiftRing;

function DOInvariantsChar2(F : PrimaryOnly := false, degmax := 10^6, degmin := 1)

    PF := Parent(F); x := PF.1; y := PF.2; z := PF.3;
    monF := [x^4, x^3*y, x^3*z, x^2*y^2, x^2*y*z, x^2*z^2, x*y^3, x*y^2*z,
	x*y*z^2, x*z^3, y^4, y^3*z, y^2*z^2, y*z^3, z^4];

    Phi := F; PolRing := BaseRing(PF); ConvFct := func<x,y|x>;
    if Characteristic(PF) ne 0 then

        PolRing := BaseRing(PF); R := LiftRing(PolRing : N := 100);

        RUVW := PolynomialRing(R, 3); U := RUVW.1; V := RUVW.2; W := RUVW.3;
        monR := [U^4, U^3*V, U^3*W, U^2*V^2, U^2*V*W, U^2*W^2, U*V^3, U*V^2*W,
            U*V*W^2, U*W^3, V^4, V^3*W, V^2*W^2, V*W^3, W^4];

        Phi := &+[ (R!ChangeBaseRing(MonomialCoefficient(F, monF[i]), R))*monR[i] : i in [1..#monF] ];

        ConvFct := ChangeBaseRing;
    end if;

    Sigma, Psi := ContravariantSigmaAndPsi(Phi);
    Rho := (1/144)*DifferentialOperation(Phi,Psi);
    He := (1/1728)*CovariantHessian(Phi);
    Tau := (1/12)*DifferentialOperation(Rho,Phi);
    Xi := (1/72)*DifferentialOperation(Sigma,He);
    Eta := (1/12)*DifferentialOperation(Xi,Sigma);
    Nu := (1/8)*DifferentialOperation(Eta,DifferentialOperation(Rho,He));

    DOs := []; WG := [];

    /* Degree 3 */
    if degmax lt 3 then return DOs, WG; end if;

    //    "I03...";
    I03 := 2^4  * 3^2 * DixmierInvariant(Phi,3 : IntegralNormalization := false);

    if degmin le 3 then
	Kx := I03;
	Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 3); /* Primary */
    end if;

    /* Degree 6 */
    if degmax lt 6 then return DOs, WG; end if;

    //    "I06...";
    I06 := 2^12 * 3^6 * DixmierInvariant(Phi,6 : IntegralNormalization := false);
    i6 := (I06 + 9*I03^2) / 8;

    if degmin le 3 then
	if not PrimaryOnly then
	    Append(~DOs, ConvFct(i6, PolRing)); Append(~WG, 6);
	end if;
    end if;

    /* Degree 9 */
    if degmax lt 9 then return DOs, WG; end if;

    //    "I09...";
    I09 := 2^12 * 3^8 * JOperation11(Tau,Rho);

    if degmin le 9 then
	Kx := I09;
	Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 9);
    end if;

    //    "J09...";
    J09 := 2^12 * 3^7 * JOperation11(Xi,Rho);

    if degmin le 9 then
	if not PrimaryOnly then
	    Kx := J09;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 9);
	end if;
    end if;

    /* Degree 12 */
    if degmax lt 12 then return DOs, WG; end if;

    //    "I12...";
    I12 := 2^16 * 3^12 * JOperation03(Rho);

    if degmin le 12 then
	Kx := I12;
	Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 12);
    end if;

    //    "J12...";
    J12 := 2^17 * 3^10 * JOperation11(Tau,Eta);
    j12 := (3*I03*I09+3*I03*J09+3*i6^2+3*I12+J12)/2;

    if degmin le 12 then
	if not PrimaryOnly then
	    Append(~DOs, ConvFct(j12, PolRing)); Append(~WG, 12);
	end if;
    end if;

    /* Degree 15 */
    if degmax lt 15 then return DOs, WG; end if;

    //    "I15...";
    I15 := 2^23 * 3^15 * JOperation30(Tau);

    if degmin le 15 then
	if not PrimaryOnly then
	    Kx := I15;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 15);
	end if;
    end if;

    //    "J15...";
    J15 := 2^23 * 3^12 * JOperation30(Xi);
    j15 :=
        (32*I03^5+32*I03^3*i6+42*I03^2*I09+60*I03^2*J09+50*I03*i6^2+26*I03*I12+36*I03*j12+48*J09*i6+43*I15+J15) / 32;

    if degmin le 15 then
	Append(~DOs, ConvFct(j15, PolRing));    Append(~WG, 15);
    end if;

    /* Degree 18 */
    if degmax lt 18 then return DOs, WG; end if;

    //    "I18...";
    I18 := 2^27 * 3^17 * JOperation22(Tau,Rho);

    if degmin le 18 then
	Kx := I18;
	Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 18);
    end if;

    //    "J18...";
    J18 := 2^27 * 3^15 * JOperation22(Xi,Rho);
    j18 :=
        (32*I03^6+24*I03^3*I09+16*I03^3*J09+40*I03^2*i6^2+16*I03^2*I12+32*I03^2*j12+8*
        I03*I09*i6+24*I03*J09*i6+16*i6^3+48*I03*I15+48*I09*J09+16*I12*i6+48*J09^2+48*
        i6*j12+47*I18+J18) / 32;

    if degmin le 18 then
	if not PrimaryOnly then
            Kx := j18;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 18);
	end if;
    end if;

    /* Degree 21 */
    if degmax lt 21 then return DOs, WG; end if;

    //    "I21...";
    I21 := 2^31 * 3^18 * JOperation03(Eta);

    if degmin le 21 then
	if not PrimaryOnly then
	    Kx := I21;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 21);
	end if;
    end if;

    //    "J21...";
    J21 := 2^33 * 3^16 * JOperation11(Nu,Eta);
    j21 :=
        (1920*I03^7+1024*I03^5*i6+684*I03^4*I09+672*I03^4*J09+1884*I03^3*i6^2+2044*I03^
        3*I12+1112*I03^3*j12+1824*I03^2*I09*i6+1088*I03^2*J09*i6+1208*I03*i6^3+389*I03
        ^2*I15+1088*I03^2*j15+8*I03*I09^2+764*I03*I09*J09+472*I03*I12*i6+572*I03*J09^2
        +1296*I03*i6*j12+176*I09*i6^2+280*J09*i6^2+1718*I03*I18+1088*I03*j18+1840*I09*
        I12+1760*I09*j12+664*I12*J09+1704*I15*i6+1520*J09*j12+1696*I21+J21) / 1024;
    if degmin le 21 then
	if not PrimaryOnly then
            Kx := j21;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 21);
	end if;
    end if;

    /* Degree 24 */
    if degmax lt 24 then return DOs, WG; end if;

    if degmin le 24 then
	if not PrimaryOnly then
	    //    "I24...";
            Kx :=
                (3*I03^2*I09*J09+3*I03^2*J09^2+3*I03*I09*i6^2+3*I03*J09*i6^2+I03^2*I18+3*I03*
                I09*I12+3*I03*I12*J09+I03*I21+3*I09*I15+3*I15*J09+3*I18*i6) / 2;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 24);
	end if;
    end if;

    if degmin le 24 then
	//    "J24...";
        Kx :=
            (6*I03^8+2*I03^6*i6+4*I03^5*I09+6*I03^5*J09+5*I03^4*i6^2+7*I03^4*I12+3*I03^4*
            j12+7*I03^3*I09*i6+3*I03^3*J09*i6+7*I03^2*i6^3+7*I03^3*I15+7*I03^3*j15+4*I03^2
            *I09^2+I03^2*I09*J09+6*I03^2*I12*i6+5*I03^2*J09^2+7*I03^2*i6*j12+7*I03*I09*i6^
            2+6*I03*J09*i6^2+7*i6^4+2*I03^2*I18+7*I03^2*j18+3*I03*I09*I12+6*I03*I09*j12+5*
            I03*I12*J09+I03*I15*i6+3*I03*J09*j12+5*I09^2*i6+7*I09*J09*i6+7*I12*i6^2+4*J09^
            2*i6+2*i6^2*j12+I03*I21+7*I03*j21+7*I09*j15+7*I12^2+2*I12*j12+4*I15*J09+6*I18*
            i6+5*J09*j15+7*j12^2) / 4;
	Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 24);
    end if;

    /* Degree 27 */
    if degmax lt 27 then return DOs, WG; end if;

    if degmin le 27 then
	//    "I27...";
	I27 := DiscriminantOfTernaryQuartic(Phi);

	Kx := I27;
	Append(~DOs, ConvFct(Kx, PolRing));     Append(~WG, 27);

	if not PrimaryOnly then

	    //    "J27...";
            Kx :=
                (6*I03^4*I15+2*I03^3*I09^2+5*I03^3*I09*J09+3*I03^3*J09^2+5*I03^2*I09*i6^2+5*I03
                ^2*J09*i6^2+7*I03^3*I18+I03^2*I09*I12+6*I03^2*I09*j12+5*I03^2*I12*J09+2*I03^2*
                I15*i6+6*I03^2*J09*j12+3*I03^2*I21+5*I03*I09*I15+3*I03*I15*J09+I03*I18*i6+2*
                I15*i6^2+2*I12*I15+6*I15*j12+6*I21*i6) / 4;
            Append(~DOs, ConvFct(Kx, PolRing));   Append(~WG, 27);

	end if;
    end if;


    if not PrimaryOnly then

	    /* Degree 30 */
	if degmax lt 30 then return DOs, WG; end if;

	if degmin le 30 then
	    //    "I30...";
            Kx :=
                (4*I03^7*I09+12*I03^6*i6^2+8*I03^5*I09*i6+12*I03^4*i6^3+14*I03^5*I15+14*I03^4*
                I09^2+5*I03^4*I09*J09+12*I03^4*I12*i6+15*I03^4*J09^2+5*I03^3*I09*i6^2+5*I03^3*
                J09*i6^2+4*I03^2*i6^4+11*I03^4*I18+5*I03^3*I09*I12+6*I03^3*I09*j12+5*I03^3*I12
                *J09+6*I03^3*I15*i6+10*I03^3*J09*j12+4*I03^2*I09^2*i6+14*I03^2*I09*J09*i6+8*
                I03^2*I12*i6^2+6*I03^2*J09^2*i6+4*I03^2*i6^2*j12+2*I03*I09*i6^3+14*I03*J09*i6^
                3+4*i6^5+7*I03^3*I21+5*I03^2*I09*I15+12*I03^2*I12^2+4*I03^2*I12*j12+3*I03^2*
                I15*J09+3*I03^2*I18*i6+2*I03*I09*I12*i6+12*I03*I09*i6*j12+14*I03*I12*J09*i6+14
                *I03*I15*i6^2+12*I03*J09*i6*j12+12*I09*J09*i6^2+12*J09^2*i6^2+12*i6^3*j12+14*
                I03*I12*I15+10*I03*I15*j12+12*I03*I21*i6+2*I09*I15*i6+4*I12^2*i6+12*I12*i6*j12
                +10*I15*J09*i6+2*I18*i6^2+8*I15^2+12*I15*j15) / 8;
	    Append(~DOs, ConvFct(Kx, PolRing));   Append(~WG, 30);

	end if;

	/* Degree 33 */
	if degmax lt 33 then return DOs, WG; end if;

	if degmin le 33 then
	    //    "I33...";
            Kx :=
                (16*I03^11+24*I03^9*i6+8*I03^8*I09+24*I03^8*J09+16*I03^7*i6^2+24*I03^7*I12+24*
                I03^7*j12+12*I03^6*I09*i6+16*I03^6*J09*i6+24*I03^5*i6^3+8*I03^6*I15+24*I03^6*
                j15+16*I03^5*I09^2+16*I03^5*I09*J09+12*I03^5*I12*i6+16*I03^5*J09^2+12*I03^5*i6
                *j12+24*I03^4*I09*i6^2+8*I03^4*J09*i6^2+8*I03^3*i6^4+12*I03^5*I18+24*I03^5*j18
                +28*I03^4*I09*I12+20*I03^4*I09*j12+12*I03^4*I12*J09+10*I03^4*I15*i6+12*I03^4*
                J09*j12+20*I03^4*i6*j15+30*I03^3*I09^2*i6+15*I03^3*I09*J09*i6+20*I03^3*I12*i6^
                2+I03^3*J09^2*i6+12*I03^3*i6^2*j12+19*I03^2*I09*i6^3+31*I03^2*J09*i6^3+8*I03*
                i6^5+24*I03^4*j21+12*I03^3*I09*I15+4*I03^3*I09*j15+16*I03^3*I12^2+16*I03^3*I12
                *j12+24*I03^3*I15*J09+17*I03^3*I18*i6+28*I03^3*J09*j15+20*I03^3*i6*j18+8*I03^3
                *j12^2+16*I03^2*I09^3+16*I03^2*I09^2*J09+19*I03^2*I09*I12*i6+28*I03^2*I09*J09^
                2+6*I03^2*I09*i6*j12+23*I03^2*I12*J09*i6+2*I03^2*I15*i6^2+28*I03^2*J09^3+30*
                I03^2*J09*i6*j12+24*I03^2*i6^2*j15+20*I03*I09^2*i6^2+20*I03*I09*J09*i6^2+24*
                I03*I12*i6^3+8*I03*J09^2*i6^2+20*I03*i6^3*j12+8*I09*i6^4+12*J09*i6^4+12*I03^2*
                I09*I18+12*I03^2*I09*j18+4*I03^2*I12*I15+24*I03^2*I12*j15+8*I03^2*I15*j12+8*
                I03^2*I18*J09+5*I03^2*I21*i6+20*I03^2*J09*j18+20*I03^2*i6*j21+7*I03*I09*I15*i6
                +12*I03*I09*J09*j12+20*I03*I09*i6*j15+20*I03*I12^2*i6+28*I03*I12*J09^2+28*I03*
                I12*i6*j12+9*I03*I15*J09*i6+7*I03*I18*i6^2+24*I03*J09^2*j12+28*I03*J09*i6*j15+
                24*I03*i6^2*j18+28*I03*i6*j12^2+28*I09^3*i6+24*I09^2*J09*i6+8*I09*I12*i6^2+16*
                I09*i6^2*j12+20*I12*J09*i6^2+30*I15*i6^3+20*J09^3*i6+16*J09*i6^2*j12+16*I03^2*
                I27+12*I03*I09*I21+20*I03*I09*j21+24*I03*I12*j18+28*I03*I15^2+28*I03*I15*j15+
                28*I03*I18*j12+4*I03*I21*J09+20*I03*J09*j21+24*I09^2*I15+20*I09^2*j15+20*I09*
                I12^2+8*I09*I12*j12+28*I09*I15*J09+20*I09*j12^2+12*I12^2*J09+22*I12*I15*i6+16*
                I12*J09*j12+16*I15*J09^2+26*I15*i6*j12+4*I18*J09*i6+22*I21*i6^2+20*J09^2*j15+
                20*J09*j12^2+20*I12*I21+28*I15*I18+20*I15*j18+28*I18*j15+20*I21*j12) / 16;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 33);
	end if;

	/* Degree 36 */
	if degmax lt 36 then return DOs, WG; end if;

	if degmin le 36 then
	    //    "I36...";
            Kx :=
                (4*I03^8*i6^2+4*I03^8*I12+4*I03^8*j12+4*I03^7*I09*i6+4*I03^7*J09*i6+2*I03^6*i6^
                3+4*I03^7*j15+4*I03^6*I09^2+2*I03^6*I09*J09+2*I03^6*I12*i6+6*I03^6*J09^2+6*I03
                ^6*i6*j12+2*I03^5*I09*i6^2+6*I03^5*J09*i6^2+2*I03^6*I18+4*I03^6*j18+4*I03^5*
                I09*I12+6*I03^5*I09*j12+4*I03^5*I12*J09+4*I03^5*I15*i6+2*I03^5*J09*j12+2*I03^5
                *i6*j15+4*I03^4*I09^2*i6+5*I03^4*I09*J09*i6+4*I03^4*I12*i6^2+I03^4*J09^2*i6+
                I03^3*I09*i6^3+3*I03^3*J09*i6^3+4*I03^2*i6^5+6*I03^5*I21+4*I03^5*j21+2*I03^4*
                I09*I15+6*I03^4*I09*j15+6*I03^4*I12^2+2*I03^4*I12*j12+6*I03^4*I15*J09+I03^4*
                I18*i6+6*I03^4*J09*j15+2*I03^4*i6*j18+4*I03^4*j12^2+6*I03^3*I09^3+I03^3*I09^2*
                J09+3*I03^3*I09*I12*i6+6*I03^3*I09*J09^2+4*I03^3*I09*i6*j12+3*I03^3*I12*J09*i6
                +3*I03^3*J09^3+2*I03^3*J09*i6*j12+2*I03^3*i6^2*j15+7*I03^2*I09^2*i6^2+2*I03^2*
                I09*J09*i6^2+I03^2*J09^2*i6^2+2*I03^2*i6^3*j12+4*I03*I09*i6^4+6*I03*J09*i6^4+
                I03^3*I09*I18+2*I03^3*I09*j18+4*I03^3*I12*I15+2*I03^3*I12*j15+6*I03^3*I15*j12+
                3*I03^3*I18*J09+3*I03^3*I21*i6+2*I03^3*J09*j18+6*I03^3*i6*j21+5*I03^2*I09^2*
                I12+6*I03^2*I09^2*j12+2*I03^2*I09*I12*J09+I03^2*I09*I15*i6+4*I03^2*I09*J09*j12
                +6*I03^2*I09*i6*j15+2*I03^2*I12^2*i6+3*I03^2*I12*J09^2+2*I03^2*I12*i6*j12+7*
                I03^2*I15*J09*i6+5*I03^2*I18*i6^2+2*I03^2*J09^2*j12+6*I03^2*J09*i6*j15+2*I03^2
                *i6^2*j18+6*I03^2*i6*j12^2+6*I03*I09^3*i6+2*I03*I09^2*J09*i6+2*I03*I09*I12*i6^
                2+4*I03*I09*J09^2*i6+2*I03*I09*i6^2*j12+2*I03*I15*i6^3+4*I03*J09*i6^2*j12+6*
                I12*i6^4+5*I03^2*I09*I21+6*I03^2*I09*j21+2*I03^2*I12*j18+4*I03^2*I15^2+2*I03^2
                *I15*j15+2*I03^2*I18*j12+3*I03^2*I21*J09+6*I03^2*J09*j21+5*I03*I09^2*I15+6*I03
                *I09^2*j15+2*I03*I09*I12*j12+5*I03*I09*I18*i6+6*I03*I09*J09*j15+6*I03*I09*j12^
                2+4*I03*I12*I15*i6+4*I03*I12*J09*j12+3*I03*I15*J09^2+4*I03*I15*i6*j12+5*I03*
                I18*J09*i6+6*I03*J09*j12^2+6*I03*i6^2*j21+6*I09^2*J09^2+2*I15*J09*i6^2+2*I18*
                i6^3+6*J09^4+6*i6^2*j12^2+6*I03*I12*I21+6*I03*I12*j21+2*I03*I15*j18+2*I03*I21*
                j12+6*I09^2*I18+4*I09*I12*I15+4*I09*I15*j12+2*I09*I18*J09+6*I12^3+4*I12*I15*
                J09+2*I12*I18*i6+6*I12*j12^2+2*I15^2*i6+2*I15*J09*j12+2*I18*J09^2+2*I18*i6*j12
                +2*I21*J09*i6+6*I15*I21+6*I15*j21+2*I18^2+6*I18*j18+6*I21*j15) / 4;

                Append(~DOs, ConvFct(Kx, PolRing));   Append(~WG, 36);
	end if;

    end if;

    return DOs, WG /*,
	[ 3, 6, 9, 9, 12, 12, 15, 15, 18, 18, 21, 21, 24, 24, 27, 27, 30, 33, 36 ] */;

    // "done...";

    return DOs, WG;

end function;
