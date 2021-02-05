/***
 *  Computes the Dixmier-Ohno invariants of a ternary quartic in char. 5
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

function DOInvariantsChar5(F : PrimaryOnly := false, degmax := 10^6, degmin := 1)

    PF := Parent(F); x := PF.1; y := PF.2; z := PF.3;
    monF := [x^4, x^3*y, x^3*z, x^2*y^2, x^2*y*z, x^2*z^2, x*y^3, x*y^2*z,
	x*y*z^2, x*z^3, y^4, y^3*z, y^2*z^2, y*z^3, z^4];

    Phi := F; PolRing := BaseRing(PF); ConvFct := func<x,y|x>;
    if Characteristic(PF) ne 0 then

        PolRing := BaseRing(PF); R := LiftRing(PolRing : N := 5);

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
    I03 := DixmierInvariant(Phi,3 : IntegralNormalization := false);

    if degmin le 3 then
	Kx := I03;
	Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 3);
    end if;

    /* Degree 6 */
    if degmax lt 6 then return DOs, WG; end if;

    //    "I06...";
    I06 := DixmierInvariant(Phi,6 : IntegralNormalization := false);

    if degmin le 6 then
	Kx := I06;
	Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 6);
    end if;

    /* Degree 9 */
    if degmax lt 9 then return DOs, WG; end if;

    //    "I09...";
    I09 := JOperation11(Tau,Rho);

    //    "J09...";
    J09 := JOperation11(Xi,Rho);

    if degmin le 9 then
	Kx := I09 + (1/5*J09 - 2/5*I09);
	Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 9);
    end if;

    if not PrimaryOnly then
	if degmin le 9 then
	    Kx := 1/5*J09 - 2/5*I09;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 9);
	end if;
    end if;

    /* Degree 12 */
    if degmax lt 12 then return DOs, WG; end if;

    //    "I12...";
    I12 := JOperation03(Rho);

    if degmin le 12 then
	Kx := I12;
	Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 12);
    end if;

    //    "J12...";
    J12 := JOperation11(Tau,Eta);

    if degmin le 12 then
	if not PrimaryOnly then
	    Kx := -1/5*I03^2*I06 - 1/5*I06^2 - 1/25*I03*J09 - 18/25*I03*I09 + 1/5*J12
		- 4/5*I12;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 12);
	end if;
    end if;

    /* Degree 15 */
    if degmax lt 15 then return DOs, WG; end if;

    //    "I15...";
    I15 := JOperation30(Tau);

    if degmin le 15 then
	if not PrimaryOnly then
	    Kx := I15;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 15);
	end if;
    end if;

    //    "J15...";
    J15 := JOperation30(Xi);

    if degmin le 15 then
	Kx := -9*I15 -64/125*I03^3*I06 - 24/125*I03*I06^2 - 79/625*I03^2*J09 - 21/125*I06*J09 +
	    253/625*I03^2*I09 + 42/125*I06*I09 - 11/125*I03*J12 + 39/125*I03*I12 +
	    1/125*J15 - 18/125*I15;
	Append(~DOs, ConvFct(Kx, PolRing));   Append(~WG, 15);
    end if;

    /* Degree 18 */
    if degmax lt 18 then return DOs, WG; end if;

    //    "I18...";
    I18 := JOperation22(Tau,Rho);

    if degmin le 18 then
	Kx := I18;
	Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 18);
    end if;

    if not PrimaryOnly then
	//    "J18...";
	J18 := JOperation22(Xi,Rho);

	if degmin le 18 then
	    Kx := 1/5*J18 - 4/5*I18;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 18);
	end if;

	/* Degree 21 */
	if degmax lt 21 then return DOs, WG; end if;

	if degmin le 21 then
	    //    "I21...";
	    I21 := JOperation03(Eta);

	    Kx := I21;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 21);
	end if;

	if degmin le 21 then
	    //    "J21...";
	    J21 := JOperation11(Nu,Eta);

	    Kx := J21;
	    Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 21);
	end if;
    end if;

    if degmin le 27 then
	//    "I27...";
	I27 := DiscriminantOfTernaryQuartic(Phi);

	Kx := I27;
	Append(~DOs, ConvFct(Kx, PolRing)); Append(~WG, 27);
    end if;

    return DOs, WG;

    // "done...";

end function;
