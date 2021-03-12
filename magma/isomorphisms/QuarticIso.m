/***
 *  Final wrapping intrinsic
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

/***
 * Exported intrinsics.
 *
 * intrinsic IsomorphismsOfTernaryQuartics(f1::RngMPolElt, f2::RngMPolElt :
 *     geometric := false) -> SeqEnum
 * intrinsic IsomorphismsOfPlaneQuartics(X1::CrvPln, X2::CrvPln :
 *     geometric := false) -> SeqEnum
 *
 * intrinsic IsIsomorphicTernaryQuartics(f1::RngMPolElt, f2::RngMPolElt :
 *     geometric := false) -> BoolElt, SeqEnum
 * intrinsic IsIsomorphicPlaneQuartics(X1::CrvPln, X2::CrvPln :
 *     geometric := false) -> BoolElt, SeqEnum
 *
 * intrinsic AutomorphismsOfTernaryQuartic(f::RngMPolElt :
 *     geometric := false) -> SeqEnum
 * intrinsic AutomorphismsOfPlaneQuartic(X::CrvPln :
 *     geometric := false) -> SeqEnum
 *
 * intrinsic AutomorphismGroupOfTernaryQuartic(f::RngMPolElt, Autos::SeqEnum :
 *     explicit := false) -> GrpPerm, Map
 * intrinsic AutomorphismGroupOfPlaneQuartic(X::CrvPln, Autos::SeqEnum :
 *     explicit := false) ->  GrpPerm, Map
 * intrinsic AutomorphismGroupOfTernaryQuartic(f::RngMPolElt :
 *     geometric := false, explicit := false) ->  GrpPerm, Map
 * intrinsic AutomorphismGroupOfPlaneQuartic(X::CrvPln :
 *     geometric := false, explicit := false) ->  GrpPerm, Map
 *
 ********************************************************************/
import "QuarticIsoFF.m": QuarticIsomorphismsFF;
import "QuarticIsoQQ.m": QuarticIsomorphismsQQ;
import "Sutherland.m": SPQIsIsomorphic;

function NormalizedM(M)

    for j := 1 to Nrows(M) do
        for i := 1 to Ncols(M) do
            if M[i,j] ne 0 then return (1 / M[i,j]) * M; end if;
        end for;
    end for;

    return M;
end function;

function ProjectiveMatrixGroup(L)

    _L := [a : a in L];
    GG, psi := GenericGroup(_L : Mult := func<a,b | NormalizedM(a*b)>);

    for i := NumberOfGenerators(GG) to 1 by -1 do
        pmp, GPrm := CosetAction(GG, sub< GG | [ GG | GG.j : j in [1..i-1] ] >);
        if #GPrm eq #GG then break; end if;
    end for;
    ReduceGenerators(~GPrm);

    return GPrm, Inverse(pmp)*psi;

end function;

/**/
intrinsic IsomorphismsOfTernaryQuartics(f1::RngMPolElt, f2::RngMPolElt :
    geometric := false) -> SeqEnum
    {Return the isomorphisms between the ternary quartic forms f1 and f2. If
    the flag geometric is set to true, then the isomorphisms over the algebraic
    closure of the base field are returned.}

    P := Parent(f1);
    require
        ((Rank(P) eq 3 and {Degree(e) : e in Monomials(f1)} eq {4}) or
        (Rank(P) eq 2 and Degree(f1) le 4))
        :
        "f1 must be a ternary homogeneous polynomial of degree 4 or a binary polynomial of degree <= 4";
    _f1 := f1;
    if Rank(P) eq 2 then
        _f1 := Basis(Homogenization(ideal<P|f1>))[1];
        if Characteristic(P) eq 0 then
            _f1 *:= LCM([Denominator(e) : e in Coefficients(_f1)]);
        end if;
    end if;

    P := Parent(f2);
    require
        ((Rank(P) eq 3 and {Degree(e) : e in Monomials(f2)} eq {4}) or
        (Rank(P) eq 2 and Degree(f2) le 4))
        :
        "f2 must be a ternary homogeneous polynomial of degree 4 or a binary polynomial of degree <= 4";
    _f2 := f2;
    if Rank(P) eq 2 then
        _f2 := Basis(Homogenization(ideal<P|f2>))[1];
        if Characteristic(P) eq 0 then
            _f2 *:= LCM([Denominator(e) : e in Coefficients(_f2)]);
        end if;
    end if;

    K := BaseRing(Parent(_f1));

    if Type(K) eq FldFin then
        _, Autos := QuarticIsomorphismsFF(_f1, _f2 : geometric := geometric);
    elif Type(K) eq FldRat then
        _, Autos := QuarticIsomorphismsQQ(_f1, _f2 : geometric := geometric);
    else
        _, Autos := SPQIsIsomorphic(_f1, _f2 : geometric := geometric);
    end if;

    return [ NormalizedM(T) : T in Autos ];

end intrinsic;

intrinsic IsomorphismsOfPlaneQuartics(X1::CrvPln, X2::CrvPln :
    geometric := false) -> SeqEnum
    {Return the isomorphisms between the plane quartic curves X1 and X2. If the
    flag geometric is set to true, then the isomorphisms over the algebraic
    closure of the base field are returned.}

    PP1 := AmbientSpace(X1);
    require IsProjective(PP1) and Dimension(PP1) eq 2 and Degree(X1) eq 4 and Genus(X1) eq 3 :
        "X1 must be a smooth projective plane quartic curve.";

    PP2 := AmbientSpace(X2);
    require IsProjective(PP2) and Dimension(PP2) eq 2 and Degree(X2) eq 4 and Genus(X2) eq 3 :
        "X2 must be a smooth projective plane quartic curve.";

    f1 := DefiningPolynomial(X1);
    f2 := DefiningPolynomial(X2);

    return IsomorphismsOfTernaryQuartics(f1, f2 : geometric := geometric);

end intrinsic;

/**/
intrinsic IsIsomorphicTernaryQuartics(f1::RngMPolElt, f2::RngMPolElt :
    geometric := false) -> BoolElt, SeqEnum
    {Determine if the ternary quartic forms f1 and f2 are isomorphic, and
    return the isomorphism between them. If the flag geometric is set to true,
    then test is performed over the algebraic closure of the base field, over
    which the isomorphisms are then determined as well.}

    Autos := IsomorphismsOfTernaryQuartics(f1, f2 : geometric := geometric);
    return #Autos ne 0, Autos;

end intrinsic;

intrinsic IsIsomorphicPlaneQuartics(X1::CrvPln, X2::CrvPln :
    geometric := false) -> BoolElt, SeqEnum
    {Determine if the plane quartic curves X1 and X2 are isomorphic, and return
    the isomorphism between them. If the flag geometric is set to true, then
    test is performed over the algebraic closure of the base field, over which
    the isomorphisms are then determined as well.}

    Autos := IsomorphismsOfPlaneQuartics(X1, X2 : geometric := geometric);
    return #Autos ne 0, Autos;

end intrinsic;

/**/
intrinsic AutomorphismsOfTernaryQuartic(f::RngMPolElt :
    geometric := false) -> SeqEnum
    {Return the automorphisms of the ternary quartic form f as matrices. If the
    flag geometric is set to true, then the automorphisms over the algebraic
    closure of the base field are returned.}

    P := Parent(f);
    require
        ((Rank(P) eq 3 and {Degree(e) : e in Monomials(f)} eq {4}) or
        (Rank(P) eq 2 and Degree(f) le 4))
        :
        "f must be a ternary homogeneous polynomial of degree 4 or a binary polynomial of degree <= 4";
    _f := f;
    if Rank(P) eq 2 then
        _f := Basis(Homogenization(ideal<P|f>))[1];
        if Characteristic(P) eq 0 then
            _f *:= LCM([Denominator(e) : e in Coefficients(_f)]);
        end if;
    end if;

    return IsomorphismsOfTernaryQuartics(_f, _f : geometric := geometric);

end intrinsic;

intrinsic AutomorphismsOfPlaneQuartic(X::CrvPln :
    geometric := false) -> SeqEnum
    {Return the automorphisms of the plane quartic curve X as matrices. If the
    flag geometric is set to true, then the automorphisms over the algebraic
    closure of the base field are returned.}

    PP := AmbientSpace(X);
    require IsProjective(PP) and Dimension(PP) eq 2 and Degree(X) eq 4 and Genus(X) eq 3 :
        "X must be a smooth projective plane quartic curve.";

    return IsomorphismsOfPlaneQuartics(X, X : geometric := geometric);

end intrinsic;

/**/
intrinsic AutomorphismGroupOfTernaryQuartic(f::RngMPolElt, Autos::SeqEnum :
    explicit := false) -> GrpPerm, Map
    {Return the automorphisms Autos of the ternary quartic form f as an
    abstract group, along with a map from said group to a matrix group if the
    flag explicit is set to true.}

    P := Parent(f);
    require
        ((Rank(P) eq 3 and {Degree(e) : e in Monomials(f)} eq {4}) or
        (Rank(P) eq 2 and Degree(f) le 4))
        :
        "f must be a ternary homogeneous polynomial of degree 4 or a binary polynomial of degree <= 4";

    aut, phi := ProjectiveMatrixGroup(Autos);

    if explicit then return aut, phi; end if;
    return aut;

end intrinsic;

intrinsic AutomorphismGroupOfPlaneQuartic(X::CrvPln, Autos::SeqEnum :
    explicit := false) ->  GrpPerm, Map
    {Return the automorphisms Autos of the plane quartic curve X as an abstract
    group, along with a map from said group to a matrix group if the flag
    explicit is set to true.}

    PP := AmbientSpace(X);
    require IsProjective(PP) and Dimension(PP) eq 2 and Degree(X) eq 4 and Genus(X) eq 3 :
        "X must be a smooth projective plane quartic curve.";

    f := DefiningPolynomial(X);
    return AutomorphismGroupOfTernaryQuartic(f, Autos : explicit := explicit);

end intrinsic;

/**/
intrinsic AutomorphismGroupOfTernaryQuartic(f::RngMPolElt :
    geometric := false, explicit := false) ->  GrpPerm, Map
    {Return the automorphisms of the ternary quartic form f as an abstract
    group, along with a map from said group to a matrix group if the flag
    explicit is set to true. If the flag geometric is set to true, then the
    automorphisms over the algebraic closure of the base field are returned.}

    _, Autos := IsIsomorphicTernaryQuartics(f, f : geometric := geometric);
    return AutomorphismGroupOfTernaryQuartic(f, Autos : explicit := explicit);

end intrinsic;

intrinsic AutomorphismGroupOfPlaneQuartic(X::CrvPln :
    geometric := false, explicit := false) ->  GrpPerm, Map
    {Return the automorphisms of the plane quartic curve X as an abstract
    group, along with a map from said group to a matrix group if the flag
    explicit is set to true. If the flag geometric is set to true, then the
    automorphisms over the algebraic closure of the base field are returned.}

    PP := AmbientSpace(X);
    require IsProjective(PP) and Dimension(PP) eq 2 and Degree(X) eq 4 and Genus(X) eq 3 :
        "X must be a smooth projective plane quartic curve.";

    f := DefiningPolynomial(X);
    return AutomorphismGroupOfTernaryQuartic(f : geometric := geometric, explicit := explicit);

end intrinsic;
