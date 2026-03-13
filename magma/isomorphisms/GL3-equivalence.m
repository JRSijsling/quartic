import "Ingredients.m":
    FindPreImage, XGCDUnique, NormalizeMatrix; 

// Return a contravariant of f using the Clebsch transfer principle
function DualTernaryForm(f) 
    d := Degree(f);
    R := Parent(f);
    Rt := PolynomialRing(R, Rank(R)-1); 

    ft := Rt!Evaluate(f, [Rt.i : i in [1..Rank(R)-1]] cat [&+[-Rt.i*R.i : i in [1..Rank(R)-1]]]);
    if d mod 2 ne 0 then
        ft := Homogenization(R!Transvectant(Transvectant(ft, ft, d-1), Transvectant(ft, ft, d-1), 2), R.Rank(R));
    else
        ft := Homogenization(R!Transvectant(ft, ft, d), R.Rank(R));
    end if;
    return ft;
end function;

/////////////////////////////////////////////////////////////////
//////////////////////// some invariants ////////////////////////
/////////////////////////////////////////////////////////////////

// Returns some invariants associated with a list of cubic covariants, with degrees given by wgt 
function SomeInvariantsTernaryCubics(L, wgt)
    return [BaseRing(Universe(L)) | TransvectantGeneral([L[1], L[2], L[i]], 3 : normalize := true) : i in [3..10]], [wgt[1]+wgt[2]+wgt[i] : i in [3..10]];
end function;


/////////////////////////////////////////////////////////////////////////
//////////////////////// defining all covariants ////////////////////////
/////////////////////////////////////////////////////////////////////////

// Return a covariant of order 4 of f, when the degree of f is not divisible by 3
function QuarticCov(f)
    assert f ne 0;
    d := Degree(f);
    assert d mod 3 ne 0;

    if d gt 15 then // decrease the bigger degrees more efficiently
        f1 := TransvectantGeneral([f, f, f], 2*Ceiling((5*d+4)/12) : normalize := true);
        d1 := Degree(f1);
        f1 := DualTernaryForm(f1);
        if d1 mod 2 eq 0 then
            f := DPairing(f1^2, f);
            d := Degree(f);
        else
            f := DPairing(f1, f);
            d := Degree(f);
        end if;
        assert f ne 0;   
    end if;
    
    if d mod 6 eq 1 then
        k := (d-1) div 6;
        c := TransvectantGeneral([f, f, f], 4*k : normalize := true);
        c1 := TransvectantGeneral([f, f, f], 4*k+2 : normalize := true);
        c2 := TransvectantGeneral([f, f, f], 6*k : normalize := true);
        return TransvectantGeneral([f, c, c1*c2], 6*k : normalize := true);
    elif d mod 6 eq 2 then
        k := (d-2) div 6;
        return TransvectantGeneral([f, f, TransvectantGeneral([f, f, f], 4*k+2 : normalize := true)], 6*k : normalize := true);
    elif d mod 6 eq 4 then
        if d eq 4 then 
            return f;
        end if;
        k := (d-4) div 6;
        c := TransvectantGeneral([f, f, f], 4*k+4 : normalize := true);
        return TransvectantGeneral([f, c, c], 6*k : normalize := true);
    elif d mod 6 eq 5 then
        k := (d-5) div 6;
        c := TransvectantGeneral([f, f, f], 4*k+4 : normalize := true);
        c1 := TransvectantGeneral([f, f, f], 6*k+4 : normalize := true);
        return TransvectantGeneral([f, f, c*c1], 6*k+4 : normalize := true);
    end if;
end function;

// Return a covariant of order 6 of f, when the degree of f is divisible by 3
function SexticCov(f)
    d := Degree(f);
    assert d mod 3 eq 0 and d ne 3;
        
    if d eq 6 then 
        return f;
    elif d mod 6 eq 0 then
        return TransvectantGeneral([f, f, f], d-2 : normalize := true);    
    end if;
    
    if d ne 9 then
        f := TransvectantGeneral([f, f, f], d-3 : normalize := true);
    end if;
    
    f1 := TransvectantGeneral([f, f, f], 8 : normalize := true);
    f1 := DualTernaryForm(f1);
    cov := DPairing(f1, f);
    cont := DPairing(cov, f1);

    return DPairing(cont, f);
end function;

function SmallOrderCovariant(f)
    d := Degree(f);
    if d mod 3 ne 0 then
        return QuarticCov(f);
    else
        return SexticCov(f);
    end if;
end function;


////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// defining the different bases of covariants ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

// Return 3 generically linearly independent covariants of order 1 of a ternary quartic 
function LinCov(f)
    assert Degree(f) eq 4;
    R := Parent(f);
    H := TransvectantGeneral([f, f, f], 2 : normalize := true); // deg 3, ord 6
    c4 := TransvectantGeneral([f, H, H], 4 : normalize := true); // deg 7, ord 4
    c5 := TransvectantGeneral([f, H, c4], 3 : normalize := true); // deg 11, ord 5
    c2 := TransvectantGeneral([f, H, c4], 4 : normalize := true); // deg 11, ord 2

    c0 := TransvectantGeneral([c5,c4,f], 4 : normalize := true); // deg 19, ord 1
    c1 := TransvectantGeneral([c5,c4,c4], 4 : normalize := true); // deg 25, ord 1
    c2 := TransvectantGeneral([c5,c2^2,f], 4 : normalize := true); // deg 34, ord 1 
    
    C := [c0,c1,c2];
    
    if Rank(Matrix([[MonomialCoefficient(c, m) : m in MonomialsOfDegree(R, 1)] : c in C])) lt #C then
        vprint Equivalence, 1 : "Warning: the covariants are not linearly independent";
        return C, [19, 25, 34], false;
    end if;
    
    return C, [19, 25, 34], true;
end function;

// Return 10 generically linearly independent covariants of order 3 of a ternary quartic 
function CubCov(f)
    assert Degree(f) eq 6;
    R := Parent(f);
    
    c206 := DualTernaryForm(f);
    
    c36 := TransvectantGeneral([f, f, f], 4 : normalize := true); // deg 3, ord 6
    c56 := TransvectantGeneral([f, f, c36], 4 : normalize := true); // deg 5, ord 6
    
    c93 := TransvectantGeneral([f, c36, c56], 5 : normalize := true); // deg 9, ord 3
    
    c1103 := DPairing(c93, c206);
    
    c123 := DPairing(c1103, f);
    c143 := DPairing(c1103, c36);
    c163 := DPairing(c1103, c56);
    
    c1403 := DPairing(c123, c206);
    
    c153 := DPairing(c1403, f);
    c173 := DPairing(c1403, c36);
    c193 := DPairing(c1403, c56);
    
    c1603 := DPairing(c143, c206);
    
    c1731 := DPairing(c1603, f);
    c1931 := DPairing(c1603, c36);
    c213 := DPairing(c1603, c56);
    
    C := [c93, c123, c143, c153, c163, c173, c1731, c193, c1931, c213];
    
    if Rank(Matrix([[MonomialCoefficient(c, m) : m in MonomialsOfDegree(R, 3)] : c in C])) lt #C then
        vprint Equivalence, 1 : "Warning: the covariants are not linearly independent";
        return C, [9, 12, 14, 15, 16, 17, 17, 19, 19, 21], false;
    end if;
    
    return C, [9, 12, 14, 15, 16, 17, 17, 19, 19, 21], true;
end function;

function BasisCovariants(f)
    d := Degree(f);
    assert d in [4, 6];
    
    if d eq 4 then
        return LinCov(f);
    end if;
    
    return CubCov(f);
end function;


///////////////////////////////////////////////////////////////
//////////////////////// main function ////////////////////////
///////////////////////////////////////////////////////////////


// Given two homogeneous polynomials in 3 variables, return true and the only possible M in GL_3 (up to scaling) such that f^M = c*g, where c is a scalar. 
// If they are not equivalent, only false is returned. When this method fails, an error is returned.
function IsGL3EquivalentFast(f, g)

    R := Parent(f);
    K := BaseRing(R);
    g := R!g;
    d := Degree(f);
    n := Rank(R);
    
    assert (n eq 3) and IsHomogeneous(f) and IsHomogeneous(g);
    assert BaseRing(Parent(f)) eq BaseRing(Parent(g));
    assert Degree(g);
    
    if d le 3 then 
        error("Not covered by this method");
    end if;    

    cov_f := SmallOrderCovariant(f);
    cov_g := SmallOrderCovariant(g);
    
    C_f, deg_basis, test_f := BasisCovariants(cov_f);
    C_g, _, test_g := BasisCovariants(cov_g);
        
    if test_f ne test_g then
        return false; 
    end if;

    if not test_f then 
        error("Not covered by this method");
    end if;
    
    if Degree(cov_f) eq 4 then
        invs_f, wgt := DixmierOhnoInvariants(cov_f);
        invs_g := DixmierOhnoInvariants(cov_g);
    else 
        invs_f, wgt := SomeInvariantsTernaryCubics(C_f, deg_basis);
        invs_g := SomeInvariantsTernaryCubics(C_g, deg_basis);
    end if;
    
    test := WPSEqual(wgt, invs_f, invs_g);

    if not test then
        return false;
    end if;
    
    if #[inv : inv in invs_f | inv eq 0] eq #invs_f then 
        error("Not covered by this method.");
    end if;

    gcd, combi := XGCDUnique([wgt[i] : i in [1..#wgt] | invs_f[i] ne 0]);
    
    if GCD([deg_basis[i]-deg_basis[i+1] : i in [1..#deg_basis-1]]) mod gcd ne 0 then
        error("Not covered by this method"); 
    end if;

    M_f := NormalizeMatrix(Transpose(Matrix(K, [[MonomialCoefficient(c, m) : m in MonomialsOfDegree(R, Degree(C_f[1]))] : c in C_f])));
    M_g := NormalizeMatrix(Transpose(Matrix(K, [[MonomialCoefficient(c, m) : m in MonomialsOfDegree(R, Degree(C_f[1]))] : c in C_g])));
       
    M_f_inv := M_f^(-1);
    
    invs_f_mod := [inv : inv in invs_f | inv ne 0];
    invs_g_mod := [inv : inv in invs_g | inv ne 0];
    inv_mult_f := &*[invs_f_mod[i]^combi[i] : i in [1..#combi]]; // normalizing factor to make all covariants of the same weight
    inv_mult_g := &*[invs_g_mod[i]^combi[i] : i in [1..#combi]];

    M_g := Transpose(Matrix([Transpose(M_g)[i]*inv_mult_g^((Max(deg_basis)-deg_basis[i]) div gcd) : i in [1..Ncols(M_g)]]));
    M_f_inv := Matrix([M_f_inv[i]/inv_mult_f^((Max(deg_basis)-deg_basis[i]) div gcd) : i in [1..Ncols(M_f)]]);

    if d mod 3 ne 0 then
        M := NormalizeMatrix(M_g*M_f_inv);
    else
        M_tilde := M_g*M_f_inv;
        M, test := FindPreImage(M_tilde, 3, 3); // now we find the actual transformation
        if not test then
            return false; // cannot be equivalent otherwise it should be defined over the base field
        end if;
        
        M := ChangeRing(Transpose(NormalizeMatrix(M)), K);
    end if;

    try 
        unit := K!(f^Transpose(M)/g);    
        return true, [* Transpose(M) *];
    catch e
        _ := K!(cov_f^Transpose(M)/cov_g); // that should always be the case, otherwise an error should be displayed
        return false; // they're not equivalent
    end try;
end function;