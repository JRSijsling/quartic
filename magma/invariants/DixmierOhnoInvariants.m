/***
 *  Computes the Dixmier-Ohno invariants of a ternary quartic.
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
 *  2016-2020 R. Lercier, C. Ritzenthaler & J.R. Sijsling
 */

 /***
 * Exported intrinsics.
 *
 * intrinsic DixmierOhnoInvariants(f::RngMPolElt : normalize := false) -> SeqEnum, SeqEnum
 *
 ********************************************************************/
import "DOmod2.m": DOInvariantsChar2;
import "DOmod3.m": DOInvariantsChar3;
import "DOmod5.m": DOInvariantsChar5;
import "DOmod7.m": DOInvariantsChar7;
import "DOmodAnyp.m": DOInvariantsCharAnyp;


function DerivativeSequence(f,n)
    S := [ Parent(f) | ];
    for k in [0..n] do
        g := f;
        for i in [1..n] do
            if i le k then
                g := Derivative(g,1);
            else
                g := Derivative(g,2);
            end if;
        end for;
        S[k+1] := g;
    end for;
    return S;
end function;


function Transvectant(f,g,r)
    P := Parent(f);
    if f eq 0 or g eq 0 then return P!0; end if;
    Sf := DerivativeSequence(f,r);
    Sg := DerivativeSequence(g,r);
    Tfg := P!0;
    for k in [0..r] do
         Tfg +:= (-1)^k*Binomial(r,k)*Sf[k+1]*Sg[r-k+1];
    end for;
    n := TotalDegree(f);
    m := TotalDegree(g);
    cfg := Factorial(n-r)*Factorial(m-r)/(Factorial(n)*Factorial(m));
    return cfg*Tfg;
end function;


function PowerDerivative(F,exp)
    DF := F;
    for i in [1..#exp] do
        if exp[i] ne 0 then
            DF := Derivative(DF,exp[i],i); // k = exp[i]-th derivative wrt i
        end if;
    end for;
    return DF;
end function;


function DifferentialOperation(F,G)
    mons := Monomials(F);
    cffs := Coefficients(F);
    DG := Parent(G)!0;
    for i in [1..#cffs] do
        DG +:= cffs[i]*PowerDerivative(G,Exponents(mons[i]));
    end for;
    return DG;
end function;


function JOperation11(F,G)
    X,Y,Z := Explode([ P.i : i in [1..3] ]) where P := Parent(F);
    ww := [  1,  1/2,  1,  1/2, 1/2,  1  ];
    XX := [ X^2, X*Y, Y^2, Y*Z, X*Z, Z^2 ];
    FF := [ MonomialCoefficient(F,M) : M in XX ];
    GG := [ MonomialCoefficient(G,M) : M in XX ];
    return &+[ ww[i]*FF[i]*GG[i] : i in [1..6] ];
end function;


function JOperation22(F,G)
    P := Parent(F);
    K := BaseRing(P);
    X,Y,Z := Explode([ P.i : i in [1..3] ]);
    A := MatrixAlgebra(K,3)!0;
    B := MatrixAlgebra(K,3)!0;
    for i in [1..3] do
        A[i,i] := MonomialCoefficient(F,P.i^2);
        B[i,i] := MonomialCoefficient(G,P.i^2);
        for j in [i+1..3] do
            A[i,j] := MonomialCoefficient(F,P.i*P.j)/2;
            A[j,i] := A[i,j];
            B[i,j] := MonomialCoefficient(G,P.i*P.j)/2;
            B[j,i] := B[i,j];
        end for;
    end for;
    Astar := Eltseq(Adjoint(A));
    Bstar := Eltseq(Adjoint(B));
    return &+[ Astar[i]*Bstar[i] : i in [1..9] ];
end function;


function JOperation30(F)
    P := Parent(F);
    K := BaseRing(P);
    X,Y,Z := Explode([ P.i : i in [1..3] ]);
    A := MatrixAlgebra(K,3)!0;
    for i in [1..3] do
        A[i,i] := MonomialCoefficient(F,P.i^2);
        for j in [i+1..3] do
            A[i,j] := MonomialCoefficient(F,P.i*P.j)/2;
            A[j,i] := A[i,j];
        end for;
    end for;
    return K!Determinant(A);
end function;


function JOperation03(G)
    return JOperation30(G);
end function;


function CovariantHessian(Phi)
    DPhi_i := [ Derivative(Phi,i) : i in [1..3] ];
    DPhi_ij := MatrixAlgebra(Parent(Phi),3)!0;
    for i in [1..3] do
        DPhi_ij[i,i] := Derivative(DPhi_i[i],i);
        for j in [i+1..3] do
            DPhi_ij[i,j] := Derivative(DPhi_i[i],j);
            DPhi_ij[j,i] := DPhi_ij[i,j];
        end for;
    end for;
    return Determinant(DPhi_ij);
end function;


function ContravariantSigmaAndPsi(Phi)

    PXYZ := Parent(Phi);
    x1, x2, x3 := Explode([ PXYZ.i : i in [1..3] ]);

    a400 := MonomialCoefficient(Phi, x1^4);
    a040 := MonomialCoefficient(Phi, x2^4);
    a004 := MonomialCoefficient(Phi, x3^4);
    a310 := MonomialCoefficient(Phi, x1^3*x2);
    a301 := MonomialCoefficient(Phi, x1^3*x3);
    a130 := MonomialCoefficient(Phi, x1*x2^3);
    a031 := MonomialCoefficient(Phi, x2^3*x3);
    a103 := MonomialCoefficient(Phi, x1*x3^3);
    a013 := MonomialCoefficient(Phi, x2*x3^3);
    a202 := MonomialCoefficient(Phi, x1^2*x3^2);
    a220 := MonomialCoefficient(Phi, x1^2*x2^2);
    a022 := MonomialCoefficient(Phi, x2^2*x3^2);
    a211 := MonomialCoefficient(Phi, x1^2*x2*x3);
    a121 := MonomialCoefficient(Phi, x1*x2^2*x3);
    a112 := MonomialCoefficient(Phi, x1*x2*x3^2);

    Sigma :=
    (a040*a004 - 1/4*a031*a013 + 1/12*a022^2)*x1^4 +
    (-a004*a130 + 1/4*a031*a103 +
    1/4*a013*a121 - 1/6*a022*a112)*x1^3*x2 + (-a040*a103 + 1/4*a130*a013 +
    1/4*a031*a112 - 1/6*a022*a121)*x1^3*x3 + (a004*a220 - 1/4*a103*a121 -
    1/4*a013*a211 + 1/6*a202*a022 + 1/12*a112^2)*x1^2*x2^2 + (3/4*a130*a103 -
    1/2*a031*a202 - 1/2*a013*a220 + 1/3*a022*a211 - 1/12*a121*a112)*x1^2*x2*x3 +
    (a040*a202 - 1/4*a130*a112 - 1/4*a031*a211 + 1/6*a220*a022 +
    1/12*a121^2)*x1^2*x3^2 + (-a004*a310 + 1/4*a301*a013 + 1/4*a103*a211 -
    1/6*a202*a112)*x1*x2^3 + (3/4*a310*a013 - 1/2*a301*a022 - 1/2*a103*a220 +
    1/3*a202*a121 - 1/12*a211*a112)*x1*x2^2*x3 + (-1/2*a310*a022 + 3/4*a301*a031
    - 1/2*a130*a202 + 1/3*a220*a112 - 1/12*a211*a121)*x1*x2*x3^2 + (-a040*a301 +
    1/4*a310*a031 + 1/4*a130*a211 - 1/6*a220*a121)*x1*x3^3 + (a400*a004 -
    1/4*a301*a103 + 1/12*a202^2)*x2^4 + (-a400*a013 + 1/4*a310*a103 +
    1/4*a301*a112 - 1/6*a202*a211)*x2^3*x3 + (a400*a022 - 1/4*a310*a112 -
    1/4*a301*a121 + 1/6*a202*a220 + 1/12*a211^2)*x2^2*x3^2 + (-a400*a031 +
    1/4*a310*a121 + 1/4*a301*a130 - 1/6*a220*a211)*x2*x3^3 + (a400*a040 -
    1/4*a310*a130 + 1/12*a220^2)*x3^4;

    Psi :=
    (1/6*a040*a004*a022 - 1/16*a040*a013^2
    - 1/16*a004*a031^2 + 1/48*a031*a013*a022
    - 1/216*a022^3)*x1^6 + (-1/6*a040*a004*a112 + 1/8*a040*a103*a013 -
    1/6*a004*a130*a022 + 1/8*a004*a031*a121 + 1/16*a130*a013^2 -
    1/48*a031*a103*a022 - 1/48*a031*a013*a112 - 1/48*a013*a022*a121 +
    1/72*a022^2*a112)*x1^5*x2 + (-1/6*a040*a004*a121 - 1/6*a040*a103*a022 +
    1/8*a040*a013*a112 + 1/8*a004*a130*a031 - 1/48*a130*a013*a022 +
    1/16*a031^2*a103 - 1/48*a031*a013*a121 - 1/48*a031*a022*a112 +
    1/72*a022^2*a121)*x1^5*x3 + (1/6*a040*a004*a202 - 1/16*a040*a103^2 +
    1/6*a004*a130*a112 - 1/8*a004*a031*a211 + 1/6*a004*a220*a022 -
    1/16*a004*a121^2 - 1/8*a130*a103*a013 + 1/48*a031*a103*a112 +
    1/48*a031*a013*a202 + 1/48*a103*a022*a121 - 1/16*a013^2*a220 +
    1/48*a013*a022*a211 + 1/48*a013*a121*a112 - 1/72*a202*a022^2 -
    1/72*a022*a112^2)*x1^4*x2^2 + (1/3*a040*a004*a211 + 1/24*a040*a103*a112 -
    1/4*a040*a013*a202 + 1/24*a004*a130*a121 - 1/4*a004*a031*a220 +
    3/16*a130*a103*a022 - 5/48*a130*a013*a112 - 5/48*a031*a103*a121 +
    1/24*a031*a013*a211 + 1/24*a031*a202*a022 + 1/48*a031*a112^2 +
    1/24*a013*a220*a022 + 1/48*a013*a121^2 - 1/36*a022^2*a211 -
    1/144*a022*a121*a112)*x1^4*x2*x3 + (1/6*a040*a004*a220 + 1/6*a040*a103*a121
    - 1/8*a040*a013*a211 + 1/6*a040*a202*a022 - 1/16*a040*a112^2 -
    1/16*a004*a130^2 - 1/8*a130*a031*a103 + 1/48*a130*a013*a121 +
    1/48*a130*a022*a112 - 1/16*a031^2*a202 + 1/48*a031*a013*a220 +
    1/48*a031*a022*a211 + 1/48*a031*a121*a112 - 1/72*a220*a022^2 -
    1/72*a022*a121^2)*x1^4*x3^2 + (-1/6*a004*a310*a022 + 1/8*a004*a301*a031 -
    1/6*a004*a130*a202 - 1/6*a004*a220*a112 + 1/8*a004*a211*a121 +
    1/16*a310*a013^2 - 1/48*a301*a013*a022 + 1/16*a130*a103^2 -
    1/48*a031*a103*a202 + 1/8*a103*a013*a220 - 1/48*a103*a022*a211 -
    1/48*a103*a121*a112 - 1/48*a013*a202*a121 - 1/48*a013*a211*a112 +
    1/36*a202*a022*a112 + 1/216*a112^3)*x1^3*x2^3 + (-1/2*a040*a004*a301 +
    1/12*a040*a103*a202 + 3/8*a004*a310*a031 - 5/24*a004*a130*a211 +
    1/12*a004*a220*a121 - 1/16*a310*a013*a022 - 1/16*a301*a031*a013 +
    1/24*a301*a022^2 - 1/16*a130*a103*a112 + 11/48*a130*a013*a202 +
    1/12*a031*a103*a211 - 1/16*a031*a202*a112 - 5/24*a103*a220*a022 +
    1/24*a103*a121^2 + 1/12*a013*a220*a112 - 1/16*a013*a211*a121 -
    1/72*a202*a022*a121 + 5/144*a022*a211*a112 - 1/144*a121*a112^2)*x1^3*x2^2*x3
    + (-1/2*a040*a004*a310 + 3/8*a040*a301*a013 - 5/24*a040*a103*a211 +
    1/12*a040*a202*a112 + 1/12*a004*a130*a220 - 1/16*a310*a031*a013 +
    1/24*a310*a022^2 - 1/16*a301*a031*a022 - 1/16*a130*a103*a121 +
    1/12*a130*a013*a211 - 5/24*a130*a202*a022 + 1/24*a130*a112^2 +
    11/48*a031*a103*a220 + 1/12*a031*a202*a121 - 1/16*a031*a211*a112 -
    1/16*a013*a220*a121 - 1/72*a220*a022*a112 + 5/144*a022*a211*a121 -
    1/144*a121^2*a112)*x1^3*x2*x3^2 + (1/8*a040*a310*a013 - 1/6*a040*a301*a022 -
    1/6*a040*a103*a220 - 1/6*a040*a202*a121 + 1/8*a040*a211*a112 -
    1/48*a310*a031*a022 + 1/16*a301*a031^2 + 1/16*a130^2*a103 +
    1/8*a130*a031*a202 - 1/48*a130*a013*a220 - 1/48*a130*a022*a211 -
    1/48*a130*a121*a112 - 1/48*a031*a220*a112 - 1/48*a031*a211*a121 +
    1/36*a220*a022*a121 + 1/216*a121^3)*x1^3*x3^3 + (1/6*a400*a004*a022 -
    1/16*a400*a013^2 + 1/6*a004*a310*a112 - 1/8*a004*a301*a121 +
    1/6*a004*a202*a220 - 1/16*a004*a211^2 - 1/8*a310*a103*a013 +
    1/48*a301*a103*a022 + 1/48*a301*a013*a112 - 1/16*a103^2*a220 +
    1/48*a103*a202*a121 + 1/48*a103*a211*a112 + 1/48*a013*a202*a211 -
    1/72*a202^2*a022 - 1/72*a202*a112^2)*x1^2*x2^4 + (-1/2*a400*a004*a031 +
    1/12*a400*a013*a022 - 5/24*a004*a310*a121 + 3/8*a004*a301*a130 +
    1/12*a004*a220*a211 + 11/48*a310*a103*a022 - 1/16*a310*a013*a112 -
    1/16*a301*a031*a103 + 1/12*a301*a013*a121 - 1/16*a301*a022*a112 -
    1/16*a130*a103*a202 + 1/24*a031*a202^2 + 1/12*a103*a220*a112 -
    1/16*a103*a211*a121 - 5/24*a013*a202*a220 + 1/24*a013*a211^2 -
    1/72*a202*a022*a211 + 5/144*a202*a121*a112 - 1/144*a211*a112^2)*x1^2*x2^3*x3
    + (a400*a040*a004 + 1/8*a400*a031*a013 - 1/12*a400*a022^2 +
    1/8*a040*a301*a103 - 1/12*a040*a202^2 + 1/8*a004*a310*a130 -
    1/12*a004*a220^2 - 5/16*a310*a031*a103 + 1/8*a310*a013*a121 -
    1/48*a310*a022*a112 - 5/16*a301*a130*a013 + 1/8*a301*a031*a112 -
    1/48*a301*a022*a121 + 1/8*a130*a103*a211 - 1/48*a130*a202*a112 -
    1/48*a031*a202*a211 - 1/48*a103*a220*a121 - 1/48*a013*a220*a211 +
    2/9*a202*a220*a022 - 5/144*a202*a121^2 - 5/144*a220*a112^2 -
    5/144*a022*a211^2 + 1/36*a211*a121*a112)*x1^2*x2^2*x3^2 +
    (-1/2*a400*a040*a013 + 1/12*a400*a031*a022 + 3/8*a040*a310*a103 -
    5/24*a040*a301*a112 + 1/12*a040*a202*a211 - 1/16*a310*a130*a013 +
    1/12*a310*a031*a112 - 1/16*a310*a022*a121 + 11/48*a301*a130*a022 -
    1/16*a301*a031*a121 - 1/16*a130*a103*a220 + 1/12*a130*a202*a121 -
    1/16*a130*a211*a112 - 5/24*a031*a202*a220 + 1/24*a031*a211^2 +
    1/24*a013*a220^2 - 1/72*a220*a022*a211 + 5/144*a220*a121*a112 -
    1/144*a211*a121^2)*x1^2*x2*x3^3 + (1/6*a400*a040*a022 - 1/16*a400*a031^2 -
    1/8*a040*a310*a112 + 1/6*a040*a301*a121 + 1/6*a040*a202*a220 -
    1/16*a040*a211^2 + 1/48*a310*a130*a022 + 1/48*a310*a031*a121 -
    1/8*a301*a130*a031 - 1/16*a130^2*a202 + 1/48*a130*a220*a112 +
    1/48*a130*a211*a121 + 1/48*a031*a220*a211 - 1/72*a220^2*a022 -
    1/72*a220*a121^2)*x1^2*x3^4 + (-1/6*a400*a004*a112 + 1/8*a400*a103*a013 -
    1/6*a004*a310*a202 + 1/8*a004*a301*a211 + 1/16*a310*a103^2 -
    1/48*a301*a103*a112 - 1/48*a301*a013*a202 - 1/48*a103*a202*a211 +
    1/72*a202^2*a112)*x1*x2^5 + (1/3*a400*a004*a121 - 1/4*a400*a103*a022 +
    1/24*a400*a013*a112 + 1/24*a004*a310*a211 - 1/4*a004*a301*a220 -
    5/48*a310*a103*a112 + 3/16*a310*a013*a202 + 1/24*a301*a103*a121 -
    5/48*a301*a013*a211 + 1/24*a301*a202*a022 + 1/48*a301*a112^2 +
    1/24*a103*a202*a220 + 1/48*a103*a211^2 - 1/36*a202^2*a121 -
    1/144*a202*a211*a112)*x1*x2^4*x3 + (-1/2*a400*a004*a130 + 3/8*a400*a031*a103
    - 5/24*a400*a013*a121 + 1/12*a400*a022*a112 + 1/12*a004*a310*a220 +
    1/12*a310*a103*a121 - 1/16*a310*a013*a211 - 5/24*a310*a202*a022 +
    1/24*a310*a112^2 - 1/16*a301*a130*a103 - 1/16*a301*a031*a202 +
    11/48*a301*a013*a220 + 1/12*a301*a022*a211 - 1/16*a301*a121*a112 +
    1/24*a130*a202^2 - 1/16*a103*a220*a211 - 1/72*a202*a220*a112 +
    5/144*a202*a211*a121 - 1/144*a211^2*a112)*x1*x2^3*x3^2 +
    (-1/2*a400*a040*a103 + 3/8*a400*a130*a013 - 5/24*a400*a031*a112 +
    1/12*a400*a022*a121 + 1/12*a040*a301*a202 - 1/16*a310*a130*a103 +
    11/48*a310*a031*a202 - 1/16*a310*a013*a220 + 1/12*a310*a022*a211 -
    1/16*a310*a121*a112 + 1/12*a301*a130*a112 - 1/16*a301*a031*a211 -
    5/24*a301*a220*a022 + 1/24*a301*a121^2 - 1/16*a130*a202*a211 +
    1/24*a103*a220^2 - 1/72*a202*a220*a121 + 5/144*a220*a211*a112 -
    1/144*a211^2*a121)*x1*x2^2*x3^3 + (1/3*a400*a040*a112 - 1/4*a400*a130*a022 +
    1/24*a400*a031*a121 - 1/4*a040*a310*a202 + 1/24*a040*a301*a211 +
    1/24*a310*a130*a112 - 5/48*a310*a031*a211 + 1/24*a310*a220*a022 +
    1/48*a310*a121^2 - 5/48*a301*a130*a121 + 3/16*a301*a031*a220 +
    1/24*a130*a202*a220 + 1/48*a130*a211^2 - 1/36*a220^2*a112 -
    1/144*a220*a211*a121)*x1*x2*x3^4 + (-1/6*a400*a040*a121 + 1/8*a400*a130*a031
    + 1/8*a040*a310*a211 - 1/6*a040*a301*a220 - 1/48*a310*a130*a121 -
    1/48*a310*a031*a220 + 1/16*a301*a130^2 - 1/48*a130*a220*a211 +
    1/72*a220^2*a121)*x1*x3^5 + (1/6*a400*a004*a202 - 1/16*a400*a103^2 -
    1/16*a004*a301^2 + 1/48*a301*a103*a202 - 1/216*a202^3)*x2^6 +
    (-1/6*a400*a004*a211 + 1/8*a400*a103*a112 - 1/6*a400*a013*a202 +
    1/8*a004*a310*a301 - 1/48*a310*a103*a202 + 1/16*a301^2*a013 -
    1/48*a301*a103*a211 - 1/48*a301*a202*a112 + 1/72*a202^2*a211)*x2^5*x3 +
    (1/6*a400*a004*a220 - 1/8*a400*a103*a121 + 1/6*a400*a013*a211 +
    1/6*a400*a202*a022 - 1/16*a400*a112^2 - 1/16*a004*a310^2 -
    1/8*a310*a301*a013 + 1/48*a310*a103*a211 + 1/48*a310*a202*a112 -
    1/16*a301^2*a022 + 1/48*a301*a103*a220 + 1/48*a301*a202*a121 +
    1/48*a301*a211*a112 - 1/72*a202^2*a220 - 1/72*a202*a211^2)*x2^4*x3^2 +
    (1/8*a400*a130*a103 - 1/6*a400*a031*a202 - 1/6*a400*a013*a220 -
    1/6*a400*a022*a211 + 1/8*a400*a121*a112 + 1/16*a310^2*a013 +
    1/8*a310*a301*a022 - 1/48*a310*a103*a220 - 1/48*a310*a202*a121 -
    1/48*a310*a211*a112 + 1/16*a301^2*a031 - 1/48*a301*a130*a202 -
    1/48*a301*a220*a112 - 1/48*a301*a211*a121 + 1/36*a202*a220*a211 +
    1/216*a211^3)*x2^3*x3^3 + (1/6*a400*a040*a202 - 1/8*a400*a130*a112 +
    1/6*a400*a031*a211 + 1/6*a400*a220*a022 - 1/16*a400*a121^2 -
    1/16*a040*a301^2 - 1/16*a310^2*a022 - 1/8*a310*a301*a031 +
    1/48*a310*a130*a202 + 1/48*a310*a220*a112 + 1/48*a310*a211*a121 +
    1/48*a301*a130*a211 + 1/48*a301*a220*a121 - 1/72*a202*a220^2 -
    1/72*a220*a211^2)*x2^2*x3^4 + (-1/6*a400*a040*a211 + 1/8*a400*a130*a121 -
    1/6*a400*a031*a220 + 1/8*a040*a310*a301 + 1/16*a310^2*a031 -
    1/48*a310*a130*a211 - 1/48*a310*a220*a121 - 1/48*a301*a130*a220 +
    1/72*a220^2*a211)*x2*x3^5 + (1/6*a400*a040*a220 - 1/16*a400*a130^2 -
    1/16*a040*a310^2 + 1/48*a310*a130*a220 - 1/216*a220^3)*x3^6;

    return Sigma, Psi;

end function;


function Rho(Phi)
    _, Psi := ContravariantSigmaAndPsi(Phi);
    return    (1/144)*DifferentialOperation(Phi, Psi);
end function;


function Tau(Phi)
    return (1/12)*DifferentialOperation(Rho(Phi), Phi);
end function;


function Xi(Phi)
    Sigma, _ := ContravariantSigmaAndPsi(Phi);
    He       := (1/1728)*CovariantHessian(Phi);
    return      (1/72)*DifferentialOperation(Sigma, He);
end function;

function DixmierInvariant(Phi,i :IntegralNormalization := false)
    P := Parent(Phi);
    K := BaseRing(P);
    if i eq 27 then
        I27 := QuarticDiscriminant(Phi);
        if IntegralNormalization then
            I27 *:= 1099511627776;
        end if;
        return I27;
    end if;
    X := P.1; Y := P.2; Z := P.3;

    a0 := MonomialCoefficient(Phi,X^4);
    b0 := (1/4)*MonomialCoefficient(Phi,X^3*Y);
    c0 := (1/6)*MonomialCoefficient(Phi,X^2*Y^2);
    d0 := (1/4)*MonomialCoefficient(Phi,Y^3*X);
    e0 := MonomialCoefficient(Phi,Y^4);
    f0 := (1/4)*MonomialCoefficient(Phi,X^3*Z);
    g0 := (1/12)*MonomialCoefficient(Phi,X^2*Y*Z);
    h0 := (1/12)*MonomialCoefficient(Phi,X*Y^2*Z);
    i0 := (1/4)*MonomialCoefficient(Phi,Y^3*Z);
    j0 := (1/6)*MonomialCoefficient(Phi,X^2*Z^2);
    k0 := (1/12)*MonomialCoefficient(Phi,X*Y*Z^2);
    l0 := (1/6)*MonomialCoefficient(Phi,Y^2*Z^2);
    m0 := (1/4)*MonomialCoefficient(Phi,X*Z^3);
    n0 := (1/4)*MonomialCoefficient(Phi,Z^3*Y);
    p0 := MonomialCoefficient(Phi,Z^4);

    if i eq 3 then
        I03 := a0*e0*p0
            + 3*(a0*l0^2+e0*j0^2+p0*c0^2)
            + 4*(b0*i0*m0+f0*d0*n0)
            - 4*(a0*i0*n0+e0*f0*m0+p0*b0*d0)
            + 6*c0*j0*l0
            + 12*(c0*k0^2+j0*h0^2+l0*g0^2)
            - 12*g0*h0*k0
            - 12*(b0*k0*l0+f0*h0*l0+d0*k0*j0+i0*g0*j0+m0*h0*c0+n0*g0*c0)
            + 12*(g0*d0*m0+h0*n0*b0+k0*f0*i0);
        if IntegralNormalization then
            I03 *:= 144;
        end if;
        return I03;
    elif i eq 6 then
        I06 := Determinant(M) where M :=Matrix(6,6, [
            a0,c0,j0,g0,f0,b0,
            c0,e0,l0,i0,h0,d0,
            j0,l0,p0,n0,m0,k0,
            g0,i0,n0,l0,k0,h0,
            f0,h0,m0,k0,j0,g0,
            b0,d0,k0,h0,g0,c0]);
        if IntegralNormalization then
            I06 *:= 2985984;
        end if;
        return I06;
    end if;

    PXYZ := Parent(Phi);
    XYZ_orig := [ PXYZ.i : i in [1..3] ];
    K := BaseRing(PXYZ);
    PUVW:= PolynomialRing(K,3); u := PUVW.1; v := PUVW.2; w := PUVW.3;
    PXYZUVW := PolynomialRing(PUVW,3); X := PXYZUVW.1; Y := PXYZUVW.2; Z := PXYZUVW.3;
    Q := Evaluate(Phi,[X,Y,Z]);
    L := u*X+v*Y+w*Z;

    R := Resultant(Q,L,Z);
    R := (-1)^Degree(Q,Z)*Resultant(Q,L,Z);

    A0 := MonomialCoefficient(R,X^4);
    B0 := (1/4)*MonomialCoefficient(R,X^3*Y);
    C0 := (1/6)*MonomialCoefficient(R,X^2*Y^2);
    D0 := (1/4)*MonomialCoefficient(R,Y^3*X);
    E0 := MonomialCoefficient(R,Y^4);

    Psi := Determinant(Matrix(3,[A0,B0,C0,B0,C0,D0,C0,D0,E0]));
    Psi := Evaluate(PUVW!Numerator(Psi/w^(TotalDegree(Psi)-6)),XYZ_orig);

    Rho := (1/144)*DifferentialOperation(Phi,Psi);
    Tau := (1/12)*DifferentialOperation(Rho,Phi);

    Xo, Yo, Zo := Explode(XYZ_orig);
    A1 := MonomialCoefficient(Tau,Xo^2);
    B1 := (1/2)*MonomialCoefficient(Tau,Xo*Yo);
    C1 := MonomialCoefficient(Tau,Yo^2);
    D1 := (1/2)*MonomialCoefficient(Tau,Xo*Zo);
    E1 := (1/2)*MonomialCoefficient(Tau,Yo*Zo);
    F1 := MonomialCoefficient(Tau,Zo^2);

    A2 := MonomialCoefficient(Rho,Xo^2);
    B2 := (1/2)*MonomialCoefficient(Rho,Xo*Yo);
    C2 := MonomialCoefficient(Rho,Yo^2);
    D2 := (1/2)*MonomialCoefficient(Rho,Xo*Zo);
    E2 := (1/2)*MonomialCoefficient(Rho,Yo*Zo);
    F2 := MonomialCoefficient(Rho,Zo^2);

    if i eq 9 then
        I09 := A1*A2+2*B1*B2+C1*C2+2*D1*D2+2*E1*E2+F1*F2;
        if IntegralNormalization then
            I09 *:= 26873856;
        end if;
        return K!I09;
    elif i eq 12 then
        I12 := Determinant(Matrix(3,3,[A2,B2,D2,B2,C2,E2,D2,E2,F2]));
        if IntegralNormalization then
            I12 *:= 34828517376;
        end if;
        return K!I12;
    elif i eq 15 then
        I15 := Determinant(Matrix(3,3,[A1,B1,D1,B1,C1,E1,D1,E1,F1]));
        if IntegralNormalization then
            I15 *:= 120367356051456;
        end if;
//      return K!MonomialCoefficient(I15,PUVW!1);
        return K!I15;
    elif i eq 18 then
        I18 := (E1^2-C1*F1)*(E2^2-C2*F2)+2*(B1*F1-D1*E1)*(B2*F2-D2*E2)
            + (D1^2-A1*F1)*(D2^2-A2*F2) + 2*(C1*D1-B1*E1)*(C2*D2-B2*E2)
            + 2*(A1*E1-B1*D1)*(A2*E2-B2*D2)+(B1^2-A1*C1)*(B2^2-A2*C2);
        if IntegralNormalization then
            I18 *:= 17332899271409664;
        end if;
//      return K!MonomialCoefficient(I18,PUVW!1);
        return K!I18;
    end if;
end function;

intrinsic CovariantHessian(Phi::RngMPolElt) -> RngMPolElt
    {}
    require Rank(Parent(Phi)) eq 3 and IsHomogeneous(Phi) :
	"Argument must be a homogeneous ternary forms.";
    DPhi_i := [ Derivative(Phi,i) : i in [1..3] ];
    DPhi_ij := MatrixAlgebra(Parent(Phi),3)!0;
    for i in [1..3] do
	DPhi_ij[i,i] := Derivative(DPhi_i[i],i);
	for j in [i+1..3] do
	    DPhi_ij[i,j] := Derivative(DPhi_i[i],j);
	    DPhi_ij[j,i] := DPhi_ij[i,j];
	end for;
    end for;
    return Determinant(DPhi_ij);
end intrinsic;


intrinsic ContravariantSigmaAndPsi(Phi::RngMPolElt) -> RngMPolElt, RngMPolElt
    {}
    // Input: Homogeneous ternary quartic.
    // Output: Contravariants Sigma and Psi of Dixmier & Ohno
    // (Salmon 3rd ed. p. 78). These should really be in the
    // dual ring PUVW, but we return them in PXYZ.
    PXYZ := Parent(Phi);
    XYZ_orig := [ PXYZ.i : i in [1..3] ];
    K := BaseRing(PXYZ);
    PUVW<u,v,w> := PolynomialRing(K,3);
    PXYZUVW<X,Y,Z> := PolynomialRing(PUVW,3);
    L := u*X + v*Y + w*Z;
    Q := Evaluate(Phi,[X,Y,Z]);
    R := Resultant(Q,L,3);
    R := (-1)^Degree(Q,Z)*Resultant(Q,L,Z);
    // This definition of Psi follows Dixmier; there is also
    // the 'symbolic notation' of Salmon, p. 78 & p. 271 as
    // \psi = (a12)^2(a23)^2(a31)^2.
    A := MonomialCoefficient(R,X^4);
    B := (1/4)*MonomialCoefficient(R,X^3*Y);
    C := (1/6)*MonomialCoefficient(R,X^2*Y^2);
    D := (1/4)*MonomialCoefficient(R,Y^3*X);
    E := MonomialCoefficient(R,Y^4);
    Psi := Determinant(Matrix(3,[A,B,C,B,C,D,C,D,E]));
    //Psi := Evaluate(Numerator(Psi/w^(TotalDegree(Psi)-6)),XYZ_orig);
    totdeg := TotalDegree(Psi)-6;
    if totdeg ge 0 then
        Psi := Evaluate(PUVW!(Psi/w^totdeg),XYZ_orig);
    else
        Psi := Evaluate(PUVW!(Psi*w^(-totdeg)),XYZ_orig);
    end if;
    // This definition of Sigma follows Salmon's 'symbolic notation' (a12)^4.
    Pxy<x,y> := PolynomialRing(PUVW,2);
    Rxy := Evaluate(R,[x,y,0]);
    Sigma := (1/2)*PUVW!Coefficients(Transvectant(Rxy,Rxy,4))[1];
    //Sigma := Evaluate(Numerator(Sigma/w^(TotalDegree(Sigma)-4)),XYZ_orig);
    totdeg := TotalDegree(Sigma)-4;
    if totdeg ge 0 then
        Sigma := Evaluate(PUVW!(Sigma/w^totdeg),XYZ_orig);
    else
        Sigma := Evaluate(PUVW!(Sigma*w^(-totdeg)),XYZ_orig);
    end if;
    return Sigma, Psi;
end intrinsic;

intrinsic QuarticCovariantsAndContravariants(Phi::RngMPolElt) -> SeqEnum
    {}
    P := Parent(Phi);
    K := BaseRing(P);
    require IsUnit(K!12) :
	"Argument must be a polynomial over a ring in which 2 and 3 are units.";
    require Rank(P) eq 3 and IsHomogeneous(Phi) and TotalDegree(Phi) eq 4 :
	"Argument must be a homogeneous ternary quartic polynomial.";
    X := P.1; Y := P.2; Z := P.3;
    // Phi itself is a contravariant:
    // Phi:   deg = 1, ord = 4
    //
    // Contravariants:
    Sigma, Psi := ContravariantSigmaAndPsi(Phi);
    Rho := (1/144)*DifferentialOperation(Phi,Psi);
    // Sigma: deg = 2, ord = 4
    // Psi:   deg = 3, ord = 6
    // Rho:   deg = 4, ord = 2
    // Covariants:
    He := (1/1728)*CovariantHessian(Phi); // deg = 3, ord = 6
    Tau := (1/12)*DifferentialOperation(Rho,Phi);
    Xi := (1/72)*DifferentialOperation(Sigma,He);
    // He:    deg = 3, ord = 6
    // Tau:   deg = 5, ord = 2
    // Xi:    deg = 5, ord = 2
    // More contravariants:
    Eta := (1/12)*DifferentialOperation(Xi,Sigma);
    Chi := (1/8)*DifferentialOperation(Tau,DifferentialOperation(Tau,Psi));
    // Eta:   deg = 7, ord = 2
    // Chi:   deg =13, ord = 2
    // More covariants:
    Nu := (1/8)*DifferentialOperation(Eta,DifferentialOperation(Rho,He));
    // Chi:   deg =14, ord = 2
    return [ Phi, He, Tau, Xi, Nu ], [ Sigma, Psi, Rho, Eta, Chi ];
end intrinsic;


intrinsic QuarticDiscriminant(f::RngMPolElt : IntegralNormalization := false) -> Any
    {Discriminant of a quartic}

    P := Parent(f);
    require
        ((Rank(P) eq 3 and {Degree(e) : e in Monomials(f)} eq {4}) or
        (Rank(P) eq 2 and Degree(f) le 4))
        :
        "Input must be a ternary homogeneous polynomial of degree  4 or a binary polynomial of degree <= 4";

    Phi := f;
    if Rank(P) eq 2 then
        Phi := Basis(Homogenization(ideal<P|f>))[1];
        if Characteristic(P) eq 0 then
            Phi *:= LCM([Denominator(e) : e in Coefficients(Phi)]);
        end if;
        P := Parent(Phi);
    end if;

    K := BaseRing(P);
    X,Y,Z := Explode([ P.i : i in [1..3] ]);

    CLASSICAL := true;

    if CLASSICAL then
        C1 := (1/4)*Derivative(Phi,X);
        C2 := (1/4)*Derivative(Phi,Y);
        C3 := (1/4)*Derivative(Phi,Z);
    else
        C1 := Derivative(Phi,X);
        C2 := Derivative(Phi,Y);
        C3 := Derivative(Phi,Z);
    end if;

    C1X2 := X^2*C1;
    C2X2 := X^2*C2;
    C3X2 := X^2*C3;
    C1Y2 := Y^2*C1;
    C2Y2 := Y^2*C2;
    C3Y2 := Y^2*C3;
    C1Z2 := Z^2*C1;
    C2Z2 := Z^2*C2;
    C3Z2 := Z^2*C3;
    C1YZ := Y*Z*C1;
    C2YZ := Y*Z*C2;
    C3YZ := Y*Z*C3;
    C1ZX := Z*X*C1;
    C2ZX := Z*X*C2;
    C3ZX := Z*X*C3;
    C1XY := X*Y*C1;
    C2XY := X*Y*C2;
    C3XY := X*Y*C3;

    if CLASSICAL then
        He := (1/1728)*CovariantHessian(Phi);
    else
        He := (1/54)*CovariantHessian(Phi);
    end if;

    if CLASSICAL then
        DHe1 := (1/2)*Derivative(He,X);
        DHe2 := (1/2)*Derivative(He,Y);
        DHe3 := (1/2)*Derivative(He,Z);
    else
        DHe1 := Derivative(He,X);
        DHe2 := Derivative(He,Y);
        DHe3 := Derivative(He,Z);
    end if;

    Eqq := [
        DHe1,DHe2,DHe3,
        C1X2,C2X2,C3X2,
        C1Y2,C2Y2,C3Y2,
        C1Z2,C2Z2,C3Z2,
        C1YZ,C2YZ,C3YZ,
        C1ZX,C2ZX,C3ZX,
        C1XY,C2XY,C3XY];

    L := [
        X^5,
        X^4*Y,
        X^4*Z,
        X^3*Y^2,
        X^3*Y*Z,
        X^3*Z^2,
        X^2*Y^3,
        X^2*Y^2*Z,
        X^2*Y*Z^2,
        X^2*Z^3,
        X*Y^4,
        X*Y^3*Z,
        X*Y^2*Z^2,
        X*Y*Z^3,
        X*Z^4,
        Y^5,
        Y^4*Z,
        Y^3*Z^2,
        Y^2*Z^3,
        Y*Z^4,
        Z^5
        ];
    R27 := Matrix(K,[ [MonomialCoefficient(Eqql,Ll): Ll in L]: Eqql in Eqq ]);
    I27 := Determinant(R27);
    if IntegralNormalization then I27 *:= 2^40; end if;

    return I27;
end intrinsic;



intrinsic DixmierOhnoInvariants(f::RngMPolElt, p::RngIntElt :
    PrimaryOnly := false,
    IntegralNormalization := false,
    degmax := Infinity(), degmin := 1,
    PolynomialOnly := false) -> SeqEnum, SeqEnum
    {
    Compute the Dixmier-Ohno invariants of a quartic given as a binary (or a ternary homogeneous) polynomial of degree 4.

    When the characteristic of the coefficient ring is 0 or greater than 7, the returned invariants are 'I3', 'I6', 'I9', 'J9', 'I12', 'J12', 'I15', 'J15', 'I18', 'J18', 'I21', 'J21' and 'I27'. If furthermore IntegralNormalization is set to true, invariants are scaled in order to be integers.

    Weights of these invariants are returned as well. If normalize is set to true, then the invariants are normalized in the corresponding weighted projective space before being returned.
    }

    P := Parent(f);
    require
        ((Rank(P) eq 3 and {Degree(e) : e in Monomials(f)} eq {4}) or
        (Rank(P) eq 2 and Degree(f) le 4))
        :
        "Input must be a ternary homogeneous polynomial of degree  4 or a binary polynomial of degree <= 4";

    F := f;
    if Rank(P) eq 2 then
        F := Basis(Homogenization(ideal<P|f>))[1];
        if Characteristic(P) eq 0 then
            F *:= LCM([Denominator(e) : e in Coefficients(F)]);
        end if;
    end if;

    if p eq 2 then
	DOs, WG := DOInvariantsChar2(F : PrimaryOnly := PrimaryOnly, degmax := degmax, degmin := degmin);

    elif p eq 3 then
	DOs, WG := DOInvariantsChar3(F : PrimaryOnly := PrimaryOnly, degmax := degmax, degmin := degmin, AllInvs := not PolynomialOnly);

    elif p eq 5 then
	DOs, WG := DOInvariantsChar5(F : PrimaryOnly := PrimaryOnly, degmax := degmax, degmin := degmin);

    elif p eq 7 then
	DOs, WG := DOInvariantsChar7(F : PrimaryOnly := PrimaryOnly, degmax := degmax, degmin := degmin);

    else
	DOs, WG := DOInvariantsCharAnyp(F : IntegralNormalization := IntegralNormalization, PrimaryOnly := PrimaryOnly, degmax := degmax, degmin := degmin);

    end if;

    return DOs, WG;

end intrinsic;

intrinsic DixmierOhnoInvariants(f::RngMPolElt :
    normalize := false,
    IntegralNormalization := false,
    PrimaryOnly := false, degmax := 10^6, degmin := 1,
    PolynomialOnly:=true) -> SeqEnum, SeqEnum
    {
    Compute the Dixmier-Ohno invariants of a quartic given as a binary (or a ternary homogeneous) polynomial of degree 4.

    When the characteristic of the coefficient ring is 0 or greater than 7, the returned invariants are 'I3', 'I6', 'I9', 'J9', 'I12', 'J12', 'I15', 'J15', 'I18', 'J18', 'I21', 'J21' and 'I27'. If furthermore IntegralNormalization is set to true, invariants are scaled in order to be integers.

    Weights of these invariants are returned as well.

    If normalize is set to true, then the invariants are normalized in the corresponding weighted projective space before being returned.
    }

    P := Parent(f);
    require
        ((Rank(P) eq 3 and {Degree(e) : e in Monomials(f)} eq {4}) or
        (Rank(P) eq 2 and Degree(f) le 4))
        :
        "Input must be a ternary homogeneous polynomial of degree  4 or a binary polynomial of degree <= 4";

    F := f;
    if Rank(P) eq 2 then
        F := Basis(Homogenization(ideal<P|f>))[1];
        if Characteristic(P) eq 0 then
            F *:= LCM([Denominator(e) : e in Coefficients(F)]);
        end if;
    end if;

    DOs, WG := DixmierOhnoInvariants(F, Characteristic(Parent(F)) :
        IntegralNormalization := IntegralNormalization,
        PrimaryOnly := PrimaryOnly,
        degmax := degmax^6, degmin := degmin,
        PolynomialOnly:= PolynomialOnly);

    if normalize eq false then return DOs, WG; end if;

    w := GCD(WG);
    return WPSNormalize([e div w : e in WG], DOs), WG;

end intrinsic;

intrinsic DixmierOhnoInvariants(C::Crv :
    normalize := false,
    IntegralNormalization := false,
    PrimaryOnly := false, degmax := 10^6, degmin := 1,
    PolynomialOnly:=true) -> SeqEnum, SeqEnum
    {
    Compute the Dixmier-Ohno invariants of a plane projective quartic curve.

    When the characteristic of the coefficient ring is 0 or greater than 7, the returned invariants are 'I3', 'I6', 'I9', 'J9', 'I12', 'J12', 'I15', 'J15', 'I18', 'J18', 'I21', 'J21' and 'I27'. If furthermore IntegralNormalization is set to true, invariants are scaled in order to be integers.

    Weights of these invariants are returned as well.

    If normalize is set to true, then the invariants are normalized in the corresponding weighted projective space before being returned.
    }

    PP := AmbientSpace(C);
    require IsProjective(PP) and Dimension(PP) eq 2 and Degree(C) eq 4 :
	"Argument must be a projective plane quartic curve.";

    return DixmierOhnoInvariants(DefiningPolynomial(C));
end intrinsic;
