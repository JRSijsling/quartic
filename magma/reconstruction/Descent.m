/***
 *  Descent functionality for plane quartics.
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
 *  2004 M. Girard, C. Ritzenthaler
 *  2006 D. Kohel
 *  2016 R. Lercier, C. Ritzenthaler & J.R. Sijsling
 */

/* TODO: We restrict to generic case for now */

import "TernaryForms.m": ConjugateForm, ConjugateMatrix;

function GL2ToGL3(U);
    a,b,c,d := Explode(Eltseq(U));
    return Matrix([[a^2, a*b, b^2], [2*a*c, a*d + b*c, 2*b*d], [c^2, c*d, d^2]]);
end function;


function IsomorphismFromB8(b8 : RandomOne := false);

    L := BaseRing(Parent(b8)); s := L.1;

    g := MinimalPolynomial(s);
    sigma := hom<L -> L | -Coefficient(g,1) - s>;

    test, L := IsGL2EquivalentExtended(b8, ConjugateForm(sigma, b8), 8 : geometric := true);

    if RandomOne then
        return GL2ToGL3(L[Random([1..#L])]);
    end if;

    return [* GL2ToGL3(ll) : ll in L *];

end function;


function NormalizeCocycle(A);

    L := BaseRing(A); s := L.1;
    g := MinimalPolynomial(s);
    sigma := hom<L -> L | -Coefficient(g,1) - s>;
    prod := A * ConjugateMatrix(sigma, A);

    lambda := [c : c in Eltseq(prod) | c ne 0][1];
    B := (lambda/Determinant(A)) * A;

    return B;

end function;


function CoboundaryRandom(A : SmallCoboundary := true, BadPrimesList := []);

    vprint QuarticRec : "";
    vprintf QuarticRec : "Looking for a coboundary with Serre's algorithm\n";

    L := BaseRing(A); s := L.1;
    g := MinimalPolynomial(s);
    sigma := hom<L -> L | -Coefficient(g,1) - s>;


    ECMLimit := 1000;
    if IsFinite(L) then
        D := L;
    else
        Bound := 10; D := [-Bound..Bound];
    end if;

    nb := 0; while true do

        // if nb gt 10 then quit; end if;

        TT := Cputime();
        nb +:= 1;
        vprint QuarticRec : "";
        vprintf QuarticRec : "\nDescent try %o :\n", nb;

        B0 := Matrix(BaseRing(A), [ [ ( Random(D) + Random(D)*s ) : c in Eltseq(row) ] : row in Rows(A) ]);
        B := B0 + A * ConjugateMatrix(sigma, B0);

        if Rank(B) eq Rank(A) then

            if not SmallCoboundary then
                return B, [];
            end if;

            nm := Norm(Determinant(B));
            ret, nm := IsSquare(nm);            /* nm is a square */

            num := Abs(Numerator(nm)); den := Abs(Denominator(nm));

            if #BadPrimesList ne 0 then
                num := num div &*[ p^Valuation(num, p) : p in BadPrimesList ];
                den := den div &*[ p^Valuation(den, p) : p in BadPrimesList ];
            end if;

            vprint QuarticRec : "";
            vprintf QuarticRec : "Checking coboundary for smallness (%o digits / %o digits)...\n", Ceiling(Log(10, num)), Ceiling(Log(10, den));

            vprint QuarticRec, 2 : "";
            vprintf QuarticRec, 2 : "den reduced := %o\n", den;

            Fac_den := Factorization(den
                : MPQSLimit := 0, ECMLimit := ECMLimit, PollardRhoLimit := 10^4, Bases := 10,  Proof := false
                );

            if (FactorizationToInteger(Fac_den) eq den) then

                vprint QuarticRec, 2 : "";
                vprintf QuarticRec, 2 : "den completely factorized, now num := %o;\n", num;

                Fac_num := Factorization(num : MPQSLimit := 0, ECMLimit
                    := ECMLimit, PollardRhoLimit := 10^4, Bases := 10, Proof := false);

                if (FactorizationToInteger(Fac_num) eq num) then
                    vprint QuarticRec, 2 : "";
                    vprintf QuarticRec, 2 : "This num is also completely factorized\n\n";

                    vprint QuarticRec : "";
                    vprintf QuarticRec : "So, we found a smooth descent morphism B, let us return it :-)\n";

                    bp := Sort([ fac[1] : fac in Fac_num ] cat [ fac[1] : fac in Fac_den ]);
                    vprint QuarticRec : "";
                    vprint QuarticRec : "Further primes in coboundary:";
                    vprint QuarticRec : bp;

                    // print "Cocycle works?";
                    //print A * ConjugateMatrix(sigma, B) eq B;

                    return B, bp;
                end if;

                vprint QuarticRec : "";
                vprintf QuarticRec :
                    "%o digits / %o digits remained... (%o s)\n\n",
                    Ceiling(Log(10, num div FactorizationToInteger(Fac_num))),
                    Ceiling(Log(10, den div
                    FactorizationToInteger(Fac_den))),
                    Cputime(TT);

            else

                vprint QuarticRec : "";
                vprintf QuarticRec :
                    "? / %o digits remained... (%o s)\n\n",
                    Ceiling(Log(10, den div FactorizationToInteger(Fac_den))),
                    Cputime(TT);

            end if;

        end if;

    end while;

end function;


function CoboundaryLinear(A : BadPrimesList := []);
    // Not great, finally...

    vprint QuarticRec : "";
    vprintf QuarticRec : "Looking for a couboundary with LLL\n";

    L := BaseRing(A); s := L.1;
    g := MinimalPolynomial(s);
    g1 := Coefficient(g,1);
    sigma := hom<L -> L | -g1 - s>;

    BB := [ Matrix(L, 3, 3, [ KroneckerDelta(i, j) : j in [1..9] ]) : i in [1..9] ]
        cat [ Matrix(L, 3, 3, [ s * KroneckerDelta(i, j) : j in [1..9] ]) : i
        in [1..9] ];

    Id := IdentityMatrix(L, 9);
    Sigma := HorizontalJoin(
        VerticalJoin(Id, -g1*Id),
        VerticalJoin(0*Id, -Id)
        );

    ASplit1 := Matrix(L, [ Eltseq(A * b) : b in BB ]);

    ASplit2 := HorizontalJoin(
        Matrix(L, [ [ Eltseq(c)[1] : c in Eltseq(row) ] : row in Rows(ASplit1)]),
        Matrix(L, [ [ Eltseq(c)[2] : c in Eltseq(row) ] : row in Rows(ASplit1) ])
        );
    K := Kernel(Sigma * ASplit2 - 1);

    if BaseRing(L) eq Rationals() then
        Lat := Lattice(Matrix(Rationals(), [ [ Rationals() ! c : c in Eltseq(K.i) ] : i in [1..Dimension(K)] ]));
        Lat := LLL(Lat);
        Lat *:= BasisDenominator(Lat);
    end if;

    Bound := 2;
    D := [-Bound..Bound];
    while true do

        v := &+[ Random(D) * Lat.i : i in [1..Rank(Lat)] ];
        B := &+[ v[i] * BB[i] : i in [1..#BB] ];
        if Rank(B) eq Rank(A) then
            //          print "Cocycle works?";
            //          print A * ConjugateMatrix(sigma, B) eq B;

            nm := Norm(Determinant(B));
            ret, nm := IsSquare(nm);        /* nm is a square */

            num := Abs(Numerator(nm));
            if #BadPrimesList ne 0 then
                num := num div &*[ p^Valuation(num, p) : p in BadPrimesList ];
            end if;

            vprint QuarticRec : "";
            vprintf QuarticRec : "Checking coboundary for smallness (%o digits)...\n", Ceiling(Log(10, num));
            vprintf QuarticRec : "num := %o\n", num;

            Fac_num := Factorization(num : MPQSLimit := 0, ECMLimit := 10^3, PollardRhoLimit := 10^4, Bases := 4, Proof := false);
            vprint QuarticRec : "";
            vprintf QuarticRec :
                "%o digits remained... \n\n",
                Ceiling(Log(10, num div FactorizationToInteger(Fac_num)));


            if (FactorizationToInteger(Fac_num) eq num) then
                vprint QuarticRec, 2 : "";
                vprintf QuarticRec, 2 : "This num is also completely factorized\n\n";

                vprint QuarticRec : "";
                vprintf QuarticRec : "So, we found a smooth descent morphism B, let us return it :-)\n";

                bp := Sort([ fac[1] : fac in Fac_num ] cat [ fac[1] : fac in Fac_num ]);
                vprint QuarticRec : "";
                vprint QuarticRec : "Further primes in coboundary:";
                vprint QuarticRec : bp;

                // print "Cocycle works?";
                //print A * ConjugateMatrix(sigma, B) eq B;

                return B, bp;
            end if;

        end if;
    end while;

end function;

function Descent(f, b8 : Isomorphism := false, RandomCoboundary := true, SmallCoboundary := false, BadPrimesList := []);

    if Type(Isomorphism) eq BoolElt then
        F := BaseRing(MinimalPolynomial(BaseRing(Parent(f)).1));

        for A in IsomorphismFromB8(b8) do
            f0 := $$(f, b8 : Isomorphism := NormalizeCocycle(A),
                RandomCoboundary := RandomCoboundary, SmallCoboundary :=
                SmallCoboundary, BadPrimesList := BadPrimesList);

            if &and( [ coeff in F : coeff in Coefficients(f0) ]) then break; end if;
        end for;

        vprint QuarticRec : "... done";
        return f0;

    end if;

    A := Isomorphism;

    if not RandomCoboundary then
        B, bp := CoboundaryLinear(A);
    else
        B, bp := CoboundaryRandom(A : SmallCoboundary := SmallCoboundary, BadPrimesList := BadPrimesList);
    end if;

    f0 := TransformForm(f, B);

    return f0 / Coefficients(f0)[1], bp;

end function;
