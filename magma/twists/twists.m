//freeze;

/***
 *  Twists of curves.
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
 *  Copyright 2020, R. Lercier, C. Ritzenthaler & J. Sijsling
 */

 /***
  *
  * Given a finite field Fq and a curve C, either hyperelliptic or a plane smooth quartic
  * Twists(C) compute a set of representatives of the twists of C. This works
  * - for C hyperelliptic g=2 and 3 without restriction on q
  * - For C hyperelliptic g>3 with q odd
  * - For C plane smooth quartic and q a power of a prime>7.
  *
  * It is based on the function TwistsOverFiniteField defined in the
  * hyperelliptic package.
  *
  *********************************************************************/

/***
 * Exported intrinsics.
 *
 * intrinsic QuarticTwists(C::Crv, Autos::SeqEnum  :
 *     AutomorphismGroup := false) -> SeqEnum[Crv], GrpPerm
 * intrinsic QuarticTwists(C::Crv :
 *     AutomorphismGroup := false) -> SeqEnum[Crv], GrpPerm
 *
 * intrinsic Twists(C::Crv :
 *     AutomorphismGroup := false) -> SeqEnum, GrpPerm
 * intrinsic GeometricAutomorphismGroup(C::Crv) -> GrpPerm
 *
 ********************************************************************/

// This function returns a Matrix with a 1 which is a multiple of M

function NormalizedM(M)

    for j := 1 to Nrows(M) do
        for i := 1 to Nrows(M) do
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


intrinsic QuarticTwists(C::Crv, Autos::SeqEnum  :
    AutomorphismGroup := false) -> SeqEnum[Crv], GrpPerm
    {Compute twisted elliptic or hyperelliptic or genus 3 plane curves, and their automorphism groups}

    F := CoefficientRing(C);

    require Type(F) eq FldFin :
        "Twist computations only available in finite fields";

    PP := AmbientSpace(C);
    require IsProjective(PP) and Dimension(PP) eq 2 and Degree(C) eq 4 and Genus(C) eq 3 :
        "C must be a smooth projective plane quartic curve.";

    Aut := [ NormalizedM(Transpose(A^(-1))) : A in Autos ];
    twists := TwistsOverFiniteField(C, Aut);
    if AutomorphismGroup then
        aut, _ := ProjectiveMatrixGroup(Aut);
        return twists, aut;
    end if;

    return twists;

end intrinsic;

intrinsic QuarticTwists(C::Crv :
    AutomorphismGroup := false) -> SeqEnum[Crv], GrpPerm
    {Compute twisted genus 3 plane curves, and their automorphism groups}

    F := CoefficientRing(C);

    require Type(F) eq FldFin :
        "Twist computations only available in finite fields";

    PP := AmbientSpace(C);
    require IsProjective(PP) and Dimension(PP) eq 2 and Degree(C) eq 4 and Genus(C) eq 3 :
        "C must be a smooth projective plane quartic curve.";

    return Twists(C : AutomorphismGroup := AutomorphismGroup);

end intrinsic;

intrinsic Twists(C::Crv :
    AutomorphismGroup := false) -> SeqEnum, GrpPerm
    {Compute twisted elliptic or hyperelliptic or genus 3 plane curves, and their geometric automorphism groups}

    F := CoefficientRing(C);

    require Type(F) eq FldFin :
	"Twist computations only available in finite fields";

    if Genus(C) eq 1 then
        Ec := EllipticCurve(C);
        twists := Twists(Ec);
        if not AutomorphismGroup then return twists; end if;
        return twists, GeometricAutomorphismGroup(Ec);
    end if;

    ishyper, H :=  IsHyperelliptic(C);
    if ishyper then
        return Twists(H : AutomorphismGroup := AutomorphismGroup);
    end if;

    PP := AmbientSpace(C);
    require IsProjective(PP) and Dimension(PP) eq 2 and Degree(C) eq 4 and Genus(C) eq 3 :
        "If not hyperelliptic, Argument must be a smooth projective plane quartic curve.";

    _, Autos := IsIsomorphicQuartic(C, C : geometric:=true);

    return QuarticTwists(C, Autos : AutomorphismGroup := AutomorphismGroup);

end intrinsic;

intrinsic GeometricAutomorphismGroup(C::Crv) -> GrpPerm
    {Compute the automorphism group of elliptic or hyperelliptic or genus 3 plane curves.}

    F := CoefficientRing(C);

    if Genus(C) eq 1 then
        Ec := EllipticCurve(C);
        return GeometricAutomorphismGroup(Ec);
    end if;

    ishyper, H :=  IsHyperelliptic(C);
    if ishyper then
        return GeometricAutomorphismGroup(H);
    end if;

    PP := AmbientSpace(C);
    require IsProjective(PP) and Dimension(PP) eq 2 and Degree(C) eq 4 and Genus(C) eq 3 :
        "If not hyperelliptic, argument must be a smooth projective plane quartic curve.";

    _, Autos := IsIsomorphicQuartic(C, C : geometric:=true);
    Autos := [ NormalizedM(Transpose(A^(-1))) : A in Autos ];
    aut, _ := ProjectiveMatrixGroup(Autos);

    return aut;

end intrinsic;
