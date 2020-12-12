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

    n := Nrows(L[1]);
    FF := BaseRing(L[1]);
    prim := PrimitiveElement(FF);
    MM := MatrixAlgebra(FF,n);
    H :=  MatrixGroup< n, FF | L cat [prim*Identity(MM)]>;
    C :=  MatrixGroup< n, FF | [prim*Identity(MM)]>;
    _, I, _ := CosetAction(H,C);

    return I;
end function;


intrinsic Twists(C::Crv : AutomorphismGroup := false) -> SeqEnum[Crv], GrpPerm
    {Compute twisted elliptic or hyperelliptic or genus 3 plane curves, and their automorphism groups}

    F := CoefficientRing(C);

    require Type(F) eq FldFin :
	"Twist computations only available in finite fields";

    if Genus(C) eq 1 then
        Ec := EllipticCurve(C);
        twists := Twists(Ec);
        if not AutomorphismGroup then return twists; end if;

        p := Characteristic(BaseRing(Ec));
        j := jInvariant(Ec);

        if j ne 0 and j ne 12^3 then
            Aut := CyclicGroup(2);
        elif p eq 2 then
            Aut := sub<Sym(24) |
                (1, 17, 9)(2, 18, 10)(3, 19, 11)(4, 20, 12)(5, 21, 13)(6, 22, 14)(7, 23, 15)(8, 24, 16),
                (1, 6, 2, 5)(3, 8, 4, 7)(9, 12, 10, 11)(13, 15, 14, 16)(17, 24, 18, 23)(19, 21, 20, 22)>;
        elif p eq 3 then
            Aut := sub<Sym(12) |
                (1, 3, 2)(4, 6, 5)(7, 8, 9)(10, 11, 12),
                (1, 10, 4, 7)(2, 11, 5, 8)(3, 12, 6, 9)>;
        elif j eq 0 then
            Aut := CyclicGroup(6);
        elif j eq 12^3 then
            Aut := CyclicGroup(4);
        end if;
        return twists, Aut;
    end if;

    ishyper, H :=  IsHyperelliptic(C);
    if ishyper then
        return Twists(H : AutomorphismGroup := AutomorphismGroup);
    end if;

    PP := AmbientSpace(C);
    require IsProjective(PP) and Dimension(PP) eq 2 and Degree(C) eq 4 and Genus(C) eq 3 :
        "Argument must be a smooth projective plane quartic curve.";

    _, Aut := IsIsomorphicQuartic(C, C : geometric:=true);
    Aut := [ NormalizedM(Transpose(A^(-1))) : A in Aut ];
    twists := TwistsOverFiniteField(C, Aut);
    if AutomorphismGroup then
        return twists, ProjectiveMatrixGroup(Aut);
    end if;

    return twists;

end intrinsic;
