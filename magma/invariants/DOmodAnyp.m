import "DixmierOhnoInvariants.m" : DerivativeSequence,  PowerDerivative, DifferentialOperation, JOperation11, JOperation22, JOperation30, JOperation03, CovariantHessian,    ContravariantSigmaAndPsi, DixmierInvariant;

function DOInvariantsCharAnyp(Phi : PrimaryOnly := false, degmax := 27, degmin := 3, p := Characteristic(BaseRing(Parent(Phi))))

    DOs := []; WG := [];

    /* Degree 3 */
    if degmax lt 3 then	return DOs, WG; end if;

    if degmin le 3 then
	//    "I03...";
	I03 := DixmierInvariant(Phi,3 : IntegralNormalization := false);

	Kx := I03;
	Append(~DOs, Kx); Append(~WG, 3);
    end if;

    /* Degree 6 */
    if degmax lt 6 then	return DOs, WG; end if;

    if degmin le 6 then
	// "I06...";
	I06 := DixmierInvariant(Phi,6 : IntegralNormalization := false);

	Kx := I06;
	Append(~DOs, Kx); Append(~WG, 6);
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
	else
	    Kx := I09;
	end if;
	Append(~DOs, Kx); Append(~WG, 9);

	if not PrimaryOnly then
	    Kx := J09;
	    Append(~DOs, Kx); Append(~WG, 9);
	end if;

    end if;

    /* Degree 12 */
    if degmax lt 12 then return DOs, WG; end if;

    Eta := (1/12)*DifferentialOperation(Xi,Sigma);

    if degmin le 12 then

	//    "I12...";
	I12 := JOperation03(Rho);

	Kx := I12;
	Append(~DOs, Kx); Append(~WG, 12);

	if not PrimaryOnly then
	    //    "J12...";
	    J12 := JOperation11(Tau,Eta);

	    Kx := J12;
	    Append(~DOs, Kx); Append(~WG, 12);
	end if;

    end if;

    /* Degree 15 */
    if degmax lt 15 then return DOs, WG; end if;

    if degmin le 15 then

	//    "I15...";
	I15 := JOperation30(Tau);
	Append(~DOs, I15); Append(~WG, 15);

	if not PrimaryOnly then
	    //    "J15...";
	    J15 := JOperation30(Xi);
	    Append(~DOs, J15); Append(~WG, 15);
	end if;

    end if;

    /* Degree 18 */
    if degmax lt 18 then return DOs, WG; end if;

    if degmin le 18 then
	//    "I18...";
	I18 := JOperation22(Tau,Rho);

	Kx := I18;
	Append(~DOs, Kx); Append(~WG, 18);


	if not PrimaryOnly then
	    //    "J18...";
	    J18 := JOperation22(Xi,Rho);

	    Kx := J18;
	    Append(~DOs, Kx); Append(~WG, 18);

	end if;

    end if;

    if not PrimaryOnly then

	/* Degree 21 */

	if degmax lt 21 then return DOs, WG; end if;

	Nu := (1/8)*DifferentialOperation(Eta,DifferentialOperation(Rho,He));

	if degmin le 21 then

	    //    "I21...";
	    I21 := JOperation03(Eta);

	    Kx := I21;
	    Append(~DOs, Kx); Append(~WG, 21);

	    //    "J21...";
	    J21 := JOperation11(Nu,Eta);

	    Kx := J21;

	    Append(~DOs, Kx); Append(~WG, 21);

	end if;

    end if;

    /* Degree 27 */

    if degmax lt 27 then return DOs, WG; end if;

    if degmin le 27 then

	//    "I27...";
	I27 := QuarticDiscriminant(Phi);

	Kx := I27;
	Append(~DOs, Kx); Append(~WG, 27);

    end if;


    return DOs, WG;

end function;
