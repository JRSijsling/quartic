/***
 *  Computes the Dixmier-Ohno invariants of a ternary quartic in char. 0 or > 7
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

function DOInvariantsCharAnyp(Phi : IntegralNormalization := false, PrimaryOnly := false, degmax := 27, degmin := 3, p := Characteristic(BaseRing(Parent(Phi))))

    DOs := []; WG := [];

    /* Degree 3 */
    if degmax lt 3 then	return DOs, WG; end if;

    if degmin le 3 then
	//    "I03...";
	I03 := DixmierInvariant(Phi,3 : IntegralNormalization := IntegralNormalization);

	Append(~DOs, I03); Append(~WG, 3);
    end if;

    /* Degree 6 */
    if degmax lt 6 then	return DOs, WG; end if;

    if degmin le 6 then
	// "I06...";
	I06 := DixmierInvariant(Phi,6 : IntegralNormalization := IntegralNormalization);

	Append(~DOs, I06); Append(~WG, 6);
    end if;

    /* Degree 9 */
    if degmax lt 9 then	return DOs, WG; end if;

    Sigma, Psi := ContravariantSigmaAndPsi(Phi);
    Rho := (1/144)*DifferentialOperation(Phi,Psi);
    Tau := (1/12)*DifferentialOperation(Rho,Phi);
    He := (1/1728)*CovariantHessian(Phi);
    Xi := (1/72)*DifferentialOperation(Sigma,He);

    if degmin le 9 then

	// "I09...";
	I09 := JOperation11(Tau,Rho);



	if not PrimaryOnly or p in {19, 47, 277, 523} then
	    //    "J09...";
	    J09 := JOperation11(Xi,Rho);
	end if;

	if p in {19, 47, 277, 523} then
	    Kx := I09 - J09;
            if IntegralNormalization then Kx := 2^12 * 3^7 * Kx; end if;
        else
            if IntegralNormalization then I09 := 2^12 * 3^8 * I09; end if;
	    Kx := I09;
	end if;
	Append(~DOs, Kx); Append(~WG, 9);

	if not PrimaryOnly then
            if IntegralNormalization then J09 := 2^12 * 3^7 * J09; end if;

	    Append(~DOs, J09); Append(~WG, 9);
	end if;

    end if;

    /* Degree 12 */
    if degmax lt 12 then return DOs, WG; end if;

    Eta := (1/12)*DifferentialOperation(Xi,Sigma);

    if degmin le 12 then

	//    "I12...";
	I12 := JOperation03(Rho);
        if IntegralNormalization then I12 := 2^16 * 3^12 * I12; end if;

	Append(~DOs, I12); Append(~WG, 12);

	if not PrimaryOnly then

            //    "J12...";
	    J12 := JOperation11(Tau,Eta);
            if IntegralNormalization then J12 := 2^17 * 3^10 * J12; end if;

	    Append(~DOs, J12); Append(~WG, 12);
	end if;

    end if;

    /* Degree 15 */
    if degmax lt 15 then return DOs, WG; end if;

    if degmin le 15 then

	//    "I15...";
	I15 := JOperation30(Tau);
        if IntegralNormalization then I15 := 2^23 * 3^15 * I15; end if;

	Append(~DOs, I15); Append(~WG, 15);

	if not PrimaryOnly then
	    //    "J15...";
	    J15 := JOperation30(Xi);
            if IntegralNormalization then J15 := 2^23 * 3^12 * J15; end if;

	    Append(~DOs, J15); Append(~WG, 15);
	end if;

    end if;

    /* Degree 18 */
    if degmax lt 18 then return DOs, WG; end if;

    if degmin le 18 then
	//    "I18...";
	I18 := JOperation22(Tau,Rho);
        if IntegralNormalization then I18 := 2^27 * 3^17 * I18; end if;

	Append(~DOs, I18); Append(~WG, 18);


	if not PrimaryOnly then
	    //    "J18...";
	    J18 := JOperation22(Xi,Rho);
            if IntegralNormalization then J18 := 2^27 * 3^15 * J18; end if;

	    Append(~DOs, J18); Append(~WG, 18);

	end if;

    end if;

    if not PrimaryOnly then

	/* Degree 21 */

	if degmax lt 21 then return DOs, WG; end if;

	Nu := (1/8)*DifferentialOperation(Eta,DifferentialOperation(Rho,He));

	if degmin le 21 then

            //    "I21...";
            I21 := JOperation03(Eta);
            if IntegralNormalization then I21 := 2^31 * 3^18 * I21; end if;

	    Append(~DOs, I21); Append(~WG, 21);

	    //    "J21...";
	    J21 := JOperation11(Nu,Eta);
            if IntegralNormalization then J21 := 2^33 * 3^16 * J21; end if;

	    Append(~DOs, J21); Append(~WG, 21);

	end if;

    end if;

    /* Degree 27 */

    if degmax lt 27 then return DOs, WG; end if;

    if degmin le 27 then

	//    "I27...";
	I27 := DiscriminantOfTernaryQuartic(Phi : IntegralNormalization := IntegralNormalization);

	Append(~DOs, I27); Append(~WG, 27);

    end if;

    //   J27 := JOperation11(Nu,Chi); // Ohno (not given name) not returned


    return DOs, WG;

end function;
