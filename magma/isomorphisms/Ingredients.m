/***
 *  Global ingredients for quartic isomorphisms
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


function FlexesThroughPoint(f0,p);

//Given a point  p  on the curve  C  defined by the smooth quartic  f0  ,
//this function determines the flex points on  C  whose tangent lines pass through  p  .

//We find these lines by simply combining the condition to be a flex line with that of passing through  p  ,
//then determining the corresponding flex points.
//This leads to overdetermination, so should be somewhat faster.

//Version zero. Perhaps also use more general algorithm with scheme structure to read off linear equation.

F := BaseRing(Parent(f0));
P2<x1,x2,x3> := ProjectiveSpace(F,2);
R2<x1,x2,x3> := CoordinateRing(P2);

f := R2!f0;

C := Curve(P2,f);
AllFlexes := Flexes(C);

if not (IsZero(F!2) or IsZero(F!3)) then


//First case:  x3 = l x1 + m x2

R<l,m,s,t> := PolynomialRing(F,4);
A4lmst := Spec(R);
A2rel<y1,y2> := AffineSpace(R,2);
P2rel<x1n,x2n,x3n> := BaseExtend(P2,R);

fn := CoordinateRing(P2rel)!f;
g := map<A2rel -> P2rel | [y1,y2,l*y1+m*y2]>;
g := AlgebraMap(g);
flm := g(fn);

a4 := R!MonomialCoefficient(flm,y1^4);
a3 := R!MonomialCoefficient(flm,y1^3*y2);
a2 := R!MonomialCoefficient(flm,y1^2*y2^2);
a1 := R!MonomialCoefficient(flm,y1*y2^3);
a0 := R!MonomialCoefficient(flm,y2^4);

//Generic flexes

p0 := p[3] - (l*p[1] + m*p[2]);
p1 := a4*s-1;
p2 := t;
p3 := 12*a4*a0 - 3*a3*a1 + a2^2;
p4 := 72*a4*a2*a0 - 27*a4*a1^2 - 27*a3^2*a0 + 9*a3*a2*a1 - 2*a2^3;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4]);
g := map<I -> P2 | [-l,-m,1]>;
Flex11 := Image(g);

//Additional set of flexes (generically empty)

p1 := s;
p2 := t;
p3 := a4;
p4 := a2^2-3*a1*a3;
p5 := a1*a2-9*a0*a3;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [-l,-m,1]>;
Flex12 := Image(g);


//Second case: x2 = l*x1

g := map<A2rel -> P2rel | [y1,l*y1,y2]>;
g := AlgebraMap(g);
flm := g(fn);

a4 := R!MonomialCoefficient(flm,y1^4);
a3 := R!MonomialCoefficient(flm,y1^3*y2);
a2 := R!MonomialCoefficient(flm,y1^2*y2^2);
a1 := R!MonomialCoefficient(flm,y1*y2^3);
a0 := R!MonomialCoefficient(flm,y2^4);

//Generic flexes

p0 := p[2] - l*p[1];
p1 := a4*s-1;
p2 := t;
p3 := 12*a4*a0 - 3*a3*a1 + a2^2;
p4 := 72*a4*a2*a0 - 27*a4*a1^2 - 27*a3^2*a0 + 9*a3*a2*a1 - 2*a2^3;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4]);
g := map<I -> P2 | [-l,1,0]>;
Flex21 := Image(g);

//Additional set of flexes (generically empty)

p1 := s;
p2 := t;
p3 := a4;
p4 := a2^2-3*a1*a3;
p5 := a1*a2-9*a0*a3;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [-l,1,0]>;
Flex22 := Image(g);


//Third case: x1 = 0

g := map<A2rel -> P2rel | [0,y1,y2]>;
g := AlgebraMap(g);
flm := g(fn);

a4 := R!MonomialCoefficient(flm,y1^4);
a3 := R!MonomialCoefficient(flm,y1^3*y2);
a2 := R!MonomialCoefficient(flm,y1^2*y2^2);
a1 := R!MonomialCoefficient(flm,y1*y2^3);
a0 := R!MonomialCoefficient(flm,y2^4);

//Generic flexes

p0 := p[1];
p1 := a4*s-1;
p2 := t;
p3 := 12*a4*a0 - 3*a3*a1 + a2^2;
p4 := 72*a4*a2*a0 - 27*a4*a1^2 - 27*a3^2*a0 + 9*a3*a2*a1 - 2*a2^3;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4]);
g := map<I -> P2 | [1,0,0]>;
Flex31 := Image(g);

//Additional set of flexes (generically empty)

p1 := s;
p2 := t;
p3 := a4;
p4 := a2^2-3*a1*a3;
p5 := a1*a2-9*a0*a3;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [1,0,0]>;
Flex32 := Image(g);


//Doing the following step, we may lose some information on Galois decompositions.

Flexdec := [Flex11,Flex12,Flex21,Flex22,Flex31,Flex32];
Flex := Flex11;
for i:=2 to #Flexdec do
	Flex := Union(Flex,Flexdec[i]);
end for;


//Characteristic 2:

elif IsZero(F!2) then


//First case:  x3 = l x1 + m x2

R<l,m,s,t> := PolynomialRing(F,4);
A4lmst := Spec(R);
A2rel<y1,y2> := AffineSpace(R,2);
P2rel<x1n,x2n,x3n> := BaseExtend(P2,R);

fn := CoordinateRing(P2rel)!f;
g := map<A2rel -> P2rel | [y1,y2,l*y1+m*y2]>;
g := AlgebraMap(g);
flm := g(fn);

a4 := R!MonomialCoefficient(flm,y1^4);
a3 := R!MonomialCoefficient(flm,y1^3*y2);
a2 := R!MonomialCoefficient(flm,y1^2*y2^2);
a1 := R!MonomialCoefficient(flm,y1*y2^3);
a0 := R!MonomialCoefficient(flm,y2^4);

// Generic flexes

p0 := p[3] - (l*p[1] + m*p[2]);
p1 := a4*s+1;
p2 := a3*t+1;
p3 := a1*a3+a2^2;
p4 := a3^2*a0+a1*(a2*a3+a1*a4);

I := Scheme(A4lmst,[p0,p1,p2,p3,p4]);
g := map<I -> P2 | [-l,-m,1]>;
Flex11 := Image(g);

//First additional set of flexes (generically empty)

p1 := s;
p2 := a4*t+1;
p3 := a3;
p4 := a2;
p5 := a1;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [-l,-m,1]>;
Flex12 := Image(g);

//Second additional set of flexes (generically empty)

p1 := s;
p2 := t;
p3 := a4;
p4 := a0*a3+a1*a2;
p5 := a1*a3+a2^2;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [-l,-m,1]>;
Flex13 := Image(g);


//Second case: x2 = l*x1

g := map<A2rel -> P2rel | [y1,l*y1,y2]>;
g := AlgebraMap(g);
flm := g(fn);

a4 := R!MonomialCoefficient(flm,y1^4);
a3 := R!MonomialCoefficient(flm,y1^3*y2);
a2 := R!MonomialCoefficient(flm,y1^2*y2^2);
a1 := R!MonomialCoefficient(flm,y1*y2^3);
a0 := R!MonomialCoefficient(flm,y2^4);

// Generic flexes

p0 := p[2] - l*p[1];
p1 := a4*s+1;
p2 := a3*t+1;
p3 := a1*a3+a2^2;
p4 := a3^2*a0+a1*(a2*a3+a1*a4);

I := Scheme(A4lmst,[p0,p1,p2,p3,p4]);
g := map<I -> P2 | [-l,1,0]>;
Flex21 := Image(g);

//First additional set of flexes (generically empty)

p1 := s;
p2 := a4*t+1;
p3 := a3;
p4 := a2;
p5 := a1;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [-l,1,0]>;
Flex22 := Image(g);

//Second additional set of flexes (generically empty)

p1 := s;
p2 := t;
p3 := a4;
p4 := a0*a3+a1*a2;
p5 := a1*a3+a2^2;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [-l,1,0]>;
Flex23 := Image(g);


//Third case: x1 = 0

g:=map<A2rel -> P2rel | [0,y1,y2]>;
g:=AlgebraMap(g);
flm:=g(fn);

a4 := R!MonomialCoefficient(flm,y1^4);
a3 := R!MonomialCoefficient(flm,y1^3*y2);
a2 := R!MonomialCoefficient(flm,y1^2*y2^2);
a1 := R!MonomialCoefficient(flm,y1*y2^3);
a0 := R!MonomialCoefficient(flm,y2^4);

// Generic flexes

p0 := p[1];
p1 := a4*s+1;
p2 := a3*t+1;
p3 := a1*a3+a2^2;
p4 := a3^2*a0+a1*(a2*a3+a1*a4);

I := Scheme(A4lmst,[p0,p1,p2,p3,p4]);
g := map<I -> P2 | [1,0,0]>;
Flex31 := Image(g);

//First additional set of flexes (generically empty)

p1 := s;
p2 := a4*t+1;
p3 := a3;
p4 := a2;
p5 := a1;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [1,0,0]>;
Flex32 := Image(g);

//Second additional set of flexes (generically empty)

p1 := s;
p2 := t;
p3 := a4;
p4 := a0*a3+a1*a2;
p5 := a1*a3+a2^2;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [1,0,0]>;
Flex33 := Image(g);


//Doing the following step, we may lose some information on Galois decompositions.

Flexdec := [Flex11,Flex12,Flex13,Flex21,Flex22,Flex23,Flex31,Flex32,Flex33];
Flex := Flex11;
for i:=2 to #Flexdec do
	Flex := Union(Flex,Flexdec[i]);
end for;


// And finally characteristic  3  .

elif IsZero(F!3) then

//First case:  x3 = l x1 + m x2

R<l,m,s,t> := PolynomialRing(F,4);
A4lmst := Spec(R);
A2rel<y1,y2> := AffineSpace(R,2);
P2rel<x1n,x2n,x3n> := BaseExtend(P2,R);

fn := CoordinateRing(P2rel)!f;
g := map<A2rel -> P2rel | [y1,y2,l*y1+m*y2]>;
g := AlgebraMap(g);
flm := g(fn);

a4 := R!MonomialCoefficient(flm,y1^4);
a3 := R!MonomialCoefficient(flm,y1^3*y2);
a2 := R!MonomialCoefficient(flm,y1^2*y2^2);
a1 := R!MonomialCoefficient(flm,y1*y2^3);
a0 := R!MonomialCoefficient(flm,y2^4);

// Generic flexes

p0 := p[3] - (l*p[1] + m*p[2]);
p1 := s;
p2 := a4*t+1;
p3 := a2;
p4 := a0*a4-a1*a3;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4]);
g := map<I -> P2 | [-l,-m,1]>;
Flex11 := Image(g);

//First additional set of flexes (generically empty)

p1 := s;
p2 := a3*t+1;
p3 := a4;
p4 := a2;
p5 := a1;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [-l,-m,1]>;
Flex12 := Image(g);

//Second additional set of flexes (generically empty)

p1 := s;
p2 := t;
p3 := a4;
p4 := a3;
p5 := a2;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [-l,-m,1]>;
Flex13 := Image(g);


//Second case: x2 = l*x1

g := map<A2rel -> P2rel | [y1,l*y1,y2]>;
g := AlgebraMap(g);
flm := g(fn);

a4 := R!MonomialCoefficient(flm,y1^4);
a3 := R!MonomialCoefficient(flm,y1^3*y2);
a2 := R!MonomialCoefficient(flm,y1^2*y2^2);
a1 := R!MonomialCoefficient(flm,y1*y2^3);
a0 := R!MonomialCoefficient(flm,y2^4);

// Generic flexes

p0 := p[2] - l*p[1];
p1 := s;
p2 := a4*t+1;
p3 := a2;
p4 := a0*a4-a1*a3;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4]);
g := map<I -> P2 | [-l,1,0]>;
Flex21 := Image(g);

//First additional set of flexes (generically empty)

p1 := s;
p2 := a3*t+1;
p3 := a4;
p4 := a2;
p5 := a1;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [-l,1,0]>;
Flex22 := Image(g);

//Second additional set of flexes (generically empty)

p1 := s;
p2 := t;
p3 := a4;
p4 := a3;
p5 := a2;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [-l,1,0]>;
Flex23 := Image(g);


//Third case: x1 = 0

g := map<A2rel -> P2rel | [0,y1,y2]>;
g := AlgebraMap(g);
flm := g(fn);

a4 := R!MonomialCoefficient(flm,y1^4);
a3 := R!MonomialCoefficient(flm,y1^3*y2);
a2 := R!MonomialCoefficient(flm,y1^2*y2^2);
a1 := R!MonomialCoefficient(flm,y1*y2^3);
a0 := R!MonomialCoefficient(flm,y2^4);

// Generic flexes

p0 := p[1];
p1 := s;
p2 := a4*t+1;
p3 := a2;
p4 := a0*a4-a1*a3;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4]);
g := map<I -> P2 | [1,0,0]>;
Flex31 := Image(g);

//First additional set of flexes (generically empty)

p1 := s;
p2 := a3*t+1;
p3 := a4;
p4 := a2;
p5 := a1;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [1,0,0]>;
Flex32 := Image(g);

//Second additional set of flexes (generically empty)

p1 := s;
p2 := t;
p3 := a4;
p4 := a3;
p5 := a2;

I := Scheme(A4lmst,[p0,p1,p2,p3,p4,p5]);
g := map<I -> P2 | [1,0,0]>;
Flex33 := Image(g);


//Doing the following step, we may lose some information on Galois decompositions.

Flexdec := [Flex11,Flex12,Flex13,Flex21,Flex22,Flex23,Flex31,Flex32,Flex33];
Flex := Flex11;
for i:=2 to #Flexdec do
	Flex := Union(Flex,Flexdec[i]);
end for;


end if;


//Test

Lns,FF := PointsOverSplittingField(Flex);
Pts := [];

P2FF<x1FF,x2FF,x3FF>:=BaseExtend(P2,FF);
fFF:=CoordinateRing(P2FF)!f;

for s in Lns do
	L := s[1]*x1FF+s[2]*x2FF+s[3]*x3FF;
	I := Scheme(P2FF,[fFF,L]);
	Ipts := PointsOverSplittingField(I);
	for pt in Ipts do
		a1 := pt[1];
		a2 := pt[2];
		a3 := pt[3];
		if (P2FF![a1,a2,a3] ne P2FF![p[1],p[2],p[3]]) then
			Append(~Pts,[a1,a2,a3]);
		end if;
	end for;
end for;

return Pts,FF;

end function;





function PowerDerivative(F,exp)
//A little copy-paste from Girard-Kohel.
    DF := F;
    for i in [1..#exp] do
	if exp[i] ne 0 then
	    DF := Derivative(DF,exp[i],i); // k = exp[i]-th derivative wrt i
	end if;
    end for;
    return DF;
end function;





function DifferentialOperation(F,G)
//A little copy-paste from Girard-Kohel.
    mons := Monomials(F);
    cffs := Coefficients(F);
    DG := Parent(G)!0;
    for i in [1..#cffs] do
	DG +:= cffs[i]*PowerDerivative(G,Exponents(mons[i]));
    end for;
    return DG;
end function;





function JOperation30(F)
//A little copy-paste from Girard-Kohel.
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





function NonSquareAsSumOfSquares(x,F)
//Expresses a non-square element of a finite field as a sum of squares.
	done := false;
	while not done do
                a := Random(F);
		y := x - a^2;
		test,b := IsSquare(y);
		if test then
			done :=true;
			soln := [a,b];
		end if;
	end while;
	return soln;
end function;





function NonSquareRepresentative(F);
//Determines a non-square in a finite field.
	done := false;
	while not done do
                a := Random(F);
		test := IsSquare(a);
		if not test then
			done := true;
		end if;
	end while;
	return a;
end function;



function RepresentativeModNthPowers(F,N);
//Determines a representative of  F^* / (F^*)^N  for a finite field  F  .
//N  is supposed to be small (since a primitive element would always work).
	M := GCD(#F-1,N);
        done := false;
	while not done do
                a := Random(F);
		test := false;
                for i:=1 to M-1 do
			test := test or IsPower(a^i,M);
		end for;
		if not test then
			done := true;
		end if;
	end while;
	return a;
end function;



function NDMSubroutine(D,n,rep,sr);
// Input:  D    diagonal matrix
//         n    single diagonal entry that is not a square
//         rep  nonsquare in  F
//         sr   precomputed square root of entry corresponding to  n  .
// Output: T    such that  T * D * T^t = Id  .
	F := Parent(D[1,1]);
	if n eq 1 then
		d2 := D[2,2];
		d3 := D[3,3];
		sr1 := sr;
		test, sr2 := IsSquare(d2*rep^(-1));
		test, sr3 := IsSquare(d3*rep^(-1));
		T1 := Matrix(F,3,3,[sr1^(-1),0,0,0,sr2^(-1),0,0,0,sr3^(-1)]);
		L := NonSquareAsSumOfSquares(rep,F);
		a := L[1];
		b := L[2];
		T2 := (1/rep)*Matrix(F,3,3,[rep,0,0,0,a,b,0,b,-a]);
		return T2 * T1;
	elif n eq 2 then
		d1 := D[1,1];
		d3 := D[3,3];
		sr2 := sr;
		test, sr1 := IsSquare(d1*rep^(-1));
		test, sr3 := IsSquare(d3*rep^(-1));
		T1 := Matrix(F,3,3,[sr1^(-1),0,0,0,sr2^(-1),0,0,0,sr3^(-1)]);
		L := NonSquareAsSumOfSquares(rep,F);
		a := L[1];
		b := L[2];
		T2 := (1/rep)*Matrix(F,3,3,[a,0,b,0,rep,0,b,0,-a]);
		return T2 * T1;
	elif n eq 3 then
		d1 := D[1,1];
		d2 := D[2,2];
		sr3 := sr;
		test, sr1 := IsSquare(d1*rep^(-1));
		test, sr2 := IsSquare(d2*rep^(-1));
		T1 := Matrix(F,3,3,[sr1^(-1),0,0,0,sr2^(-1),0,0,0,sr3^(-1)]);
		L := NonSquareAsSumOfSquares(rep,F);
		a := L[1];
		b := L[2];
		T2 := (1/rep)*Matrix(F,3,3,[a,b,0,b,-a,0,0,0,rep]);
		return T2 * T1;
	end if;
end function;





function NormalizeDiagonalMatrix(D,rep)
// Input:  D    diagonal matrix
//         rep  nonsquare in F
// Output: T    such that either  T * D * T^t = Id  or  T * D * T^t = rep * Id  .
	F := Parent(D[1,1]);
	d1 := D[1,1];
	d2 := D[2,2];
	d3 := D[3,3];
	nSq := 0;
	is1,sr1 := IsSquare(d1);
	if is1 then
		nSq := nSq + 1;
	end if;
	is2,sr2 := IsSquare(d2);
	if is2 then
		nSq := nSq + 1;
	end if;
	is3,sr3 := IsSquare(d3);
	if is3 then
		nSq := nSq + 1;
	end if;
	if nSq eq 3 then
		T := Matrix(F,3,3,[sr1^(-1),0,0,0,sr2^(-1),0,0,0,sr3^(-1)]);
		return [IdentityMatrix(F,3) , T];
	elif nSq eq 2 then
		T := NormalizeDiagonalMatrix(rep^(-1)*D,rep)[2];
		return [rep * IdentityMatrix(F,3) , T];
	elif nSq eq 1 then
		if is1 then
			return [IdentityMatrix(F,3), NDMSubroutine(D,1,rep,sr1)];
		elif is2 then
			return [IdentityMatrix(F,3), NDMSubroutine(D,2,rep,sr2)];
		elif is3 then
			return [IdentityMatrix(F,3), NDMSubroutine(D,3,rep,sr3)];
		end if;
	elif nSq eq 0 then
		T := NormalizeDiagonalMatrix(rep^(-1)*D,rep)[2];
		return [rep * IdentityMatrix(F,3) , T];
	end if;
end function;





function ThirdRootModFourthPowers(x,n);
// Third root mod 4n-th powers.

a,b,c := ExtendedGreatestCommonDivisor(3,4*n);

return x^b;

end function;





function SmallSplittingFieldOverRationals(D);
//Given discriminant of a quaternion algebra, determines a small splitting field

if D eq 1 then
	return 1;
else

Fac := Factorization(D);
d := -1;
done := false;
while not done do
	if IsSquarefree(d) then
		done := true;
		if IsOdd(#Fac) then
			d := -d;
		end if;
		for fac in Fac do
			if IsEven(fac[1]) then
				if (d mod 8) eq 1 then
					done :=false;
				end if;
			else
				if LegendreSymbol(d,fac[1]) eq 1 then
					done:=false;
				end if;
			end if;
		end for;
	end if;
	if not done then
		d := AbsoluteValue(d);
		d := d + 1;
	end if;
end while;

return d;

end if;

end function;





function EffSPProduct(a,b)
//Efficiently determines the square part of a product of two rational numbers.
//Minus is neglected here, but that does not matter for squarefree part.

FacNa := Factorization(Numerator(a));
FacDa := Factorization(Denominator(a));
FacNb := Factorization(Numerator(b));
FacDb := Factorization(Denominator(b));

FacN := FacNa*FacNb;
FacD := FacDa*FacDb;

GCDND := GCD(FacN,FacD);
FacN := FacN/GCDND;
FacD := FacD/GCDND;

xN,yN := SquareFreeFactorization(FacN);
xD,yD := SquareFreeFactorization(FacD);

yN := Integers()!yN;
xD := Integers()!xD;
yD := Integers()!yD;

return yN/(xD*yD);

end function;





function EffSPProduct(a,b)
//Same as above, but avoids complete factorization.
//May miss canceling factors.

Na := Numerator(a);
Da := Denominator(a);
Nb := Numerator(b);
Db := Denominator(b);

xNa,yNa := SquareFreeFactorization(Na);
xNb,yNb := SquareFreeFactorization(Nb);
xDa,yDa := SquareFreeFactorization(Da);
xDb,yDb := SquareFreeFactorization(Db);
xPN,yPN := SquareFreeFactorization(xNa*xNb);
xPD,yPD := SquareFreeFactorization(xDa*xDb);

yN := Integers()!yPN*yNa*yNb;
xD := Integers()!xPD;
yD := Integers()!yPD*yDa*yDb;

return yN/(xD*yD);

end function;





function TupleToSequence(tup);

L:=[];
for i:=1 to #tup do
	Append(~L,tup[i]);
end for;

return L;

end function;





function SmallSplittingFieldNaive(Dfin,Dinf,F);
//Given discriminant  D  of a quaternion algebra over a totally real field  F
//(as an ideal),
//determines a small splitting field for the algebra.
//This field is return in the form of an integral element
//of which a square root is to be drawn.

if (#Dfin eq 1) and (#Dinf eq 0) then
	return F;
else

ZF := Integers(F);
deg := Degree(F);

lim := 1;
found := false;

while not found do
	Indices := CartesianPower([-lim..lim],deg);
	for index in Indices do
		if not found then
			seq := TupleToSequence(index);
			d := ZF!seq;
			if d ne 0 then
				test := true;
				for dinf in Dinf do
					if test then
						test := Evaluate(d,dinf) lt 0;
					end if;
				end for;
				for dfin in Dfin do
					if test then
						rcf,rcfm := ResidueClassField(dfin);
						res := rcfm(d);
						//if IsEven(Norm(dfin)) then
						//	test := IsOdd(Valuation(d,dfin));
						//else
							test := (IsOdd(Valuation(d,dfin)) or (not IsSquare(res)));
						//end if;
					end if;
				end for;
				if test then
					found := true;
				end if;
			end if;
		end if;
	end for;
	lim := lim + 1;
end while;

R<t> := PolynomialRing(F);
return NumberField(t^2-d);

end if;

end function;





function TrivializeAlgebra(B,pt)
//Given a quaternion algebra  B  with trivial discriminant,
//and an element of norm zero,
//this function determines a trivialization of  B  .
//The ideas used are due to John Voight.

F := BaseField(B);
//A := AbsoluteField(F);

i := B![0,1,0,0];
j := B![0,0,1,0];
k := B![0,0,0,1];

//a := A!(F!(-i^2));
//b := A!(F!(-j^2));
//c := A!(F!(-k^2));

//P2<x,y,z> := ProjectiveSpace(A,2);
//C := Conic(P2,a*x^2 + b*y^2 + c*z^2);

//test,pt := HasRationalPoint(C);
//pt := [F!pt[1],F!pt[2],F!pt[3]];

e := pt[1]*i + pt[2]*j + pt[3]*k;

l := i;
s := Trace(e*i);
if s eq 0 then
	s := Trace(e*j);
	if s eq 0 then
		s := Trace(e*k);
		l := k;
	else
		l := j;
	end if;
end if;

t := Trace(l);
n := Norm(l);
ep := -(1/s)*e;

jp := l + (t*l - n - 1)*ep;
if Characteristic(F) ne 2 then
	ip := -(ep*l - (l - t)*ep);
else
	ip := l + ((t + 1)*l + n + 1)*ep;
end if;

return ip,jp;

end function;





function Normalize22(M);
//Divides out scalars in 2x2 invertible matrix.

a := M[1,1];
b := M[1,2];
c := M[2,1];
d := M[2,2];
F := Parent(a);

if a ne 0 then
	return Matrix(F,2,2,[1,b/a,c/a,d/a]);
else
	return Matrix(F,2,2,[0,b/c,1,d/c]);
end if;

end function;





function Normalize33(M);
//Divides out scalars in 3x3 invertible matrix.

S := Eltseq(M);
i0 := Minimum([ i : i in [1..#S] | S[i] ne 0 ]);
return (1/S[i0])*M;

end function;





function IsMultiple(f,g);
//Determines if the form  f  is a multiple of the non-zero form  g  .

Cf := Coefficients(f);
Cg := Coefficients(g);

stop := false;
for i := 1 to #Cg do
	if not stop then
		if Cg[i] ne 0 then
			factor := Cf[i]/Cg[i];
			stop := true;
		end if;
	end if;
end for;

test := f eq factor*g;

return test,factor;

end function;





function ConjugatePol(f,sigma);

Mons := Monomials(f);
fc := 0;

for mon in Mons do
    fc := fc + sigma(MonomialCoefficient(f,mon))*mon;
end for;

return fc;

end function;





function ConjugatePolFrob(f,r);

Mons := Monomials(f);
fc := 0;

for mon in Mons do
    fc := fc + Frobenius(MonomialCoefficient(f,mon),r)*mon;
end for;

return fc;

end function;





function ConjugateMat(M,sigma);

MRing := Parent(M);
CM := ElementToSequence(M);
CMc := [];

for i:=1 to #CM do
    CMc[i] := sigma(CM[i]);
end for;

return MRing!CMc;

end function;





function ConjugateMatFrob(M,r);

MRing := Parent(M);
CM := ElementToSequence(M);
CMc := [];

for i:=1 to #CM do
    CMc[i] := Frobenius(CM[i],r);
end for;

return MRing!CMc;

end function;





function BinQuadInvs(bq);

R<x,z> := Parent(bq);

a := MonomialCoefficient(bq,x^4);
b := MonomialCoefficient(bq,x^3*z);
c := MonomialCoefficient(bq,x^2*z^2);
d := MonomialCoefficient(bq,x*z^3);
e := MonomialCoefficient(bq,z^4);

I := 12*a*e - 3*b*d + c^2;
J := 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3;
Delta := 4*I^3 - J^2;

return I,J,Delta;

end function;





function TransformBinaryForm(f,T);
//Input :   A homogeneous binary form  f  and an invertible matrix  T  over the same ground field.
//Output :  The transformation of  f  by  T  . Contravariant.

S2 := Parent(f);
x := S2.1;
z := S2.2;

//Ti := T^(-1);
Ti := T;
h := hom<S2 -> S2 | [Ti[1,1]*x + Ti[1,2]*z, Ti[2,1]*x + Ti[2,2]*z]>;

return h(f);

end function;





function TransformTernaryForm(f,T);
//Input :   A homogeneous ternary form  f  and an invertible matrix  T  over the same ground field.
//Output :  The transformation of  f  by  T  . Contravariant.

S3 := Parent(f);
x := S3.1;
y := S3.2;
z := S3.3;

//Ti := T^(-1);
Ti := T;
h := hom<S3 -> S3 | [Ti[1,1]*x + Ti[1,2]*y + Ti[1,3]*z, Ti[2,1]*x + Ti[2,2]*y + Ti[2,3]*z, Ti[3,1]*x + Ti[3,2]*y + Ti[3,3]*z]>;

return h(f);

end function;





function IsMultiplePolynomial(f, f0)
// Returns whether f is a multiple of f0 or not, along with a corresponding scalar if this is the case.

if not Monomials(f) eq Monomials(f0) then
    return false, 0;
end if;
coeffs := Coefficients(f); coeffs0 := Coefficients(f0);
M := Matrix([ coeffs ]); M0 := Matrix([ coeffs0 ]);
test, A := IsConsistent(M0, M);
if not test then
    return false, 0;
else
    return true, A[1,1];
end if;

end function;





procedure AssertTs(f1, f2, Ts : geometric := false)

if not geometric then
    for T in Ts do
        assert IsMultiplePolynomial(TransformForm(f1,T), f2);
    end for;
else
    for T in Ts do
        R := Parent(f1); F := BaseRing(R);
        K := BaseRing(T); S := PolynomialRing(K, 3);
        /* We rely on coercion here */
        h := hom< R -> S | [S.1,S.2,S.3] >;
        g1 := h(f1); g2 := h(f2);
        assert IsMultiplePolynomial(TransformForm(g1,T), g2);
    end for;
end if;

end procedure;
