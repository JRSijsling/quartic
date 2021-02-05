// Non-hyperelliptic case

NormalizedM := function(M)
  N := Nrows(M);
  for j in [1 .. N] do
    for i in [1 .. N] do
      if M[i,j] ne 0 then
        return 1/M[i,j]*M;
      end if;
    end for;
  end for;
  printf("Error, your matrix is zero");
  return 0;
end function;

P<x,y,z>:=PolynomialRing(GF(5),3);
PP:=ProjectiveSpace(P);
f:=x^4+y^4+z^4+3*(x^2*y^2+y^2*z^2);
//f:=x^3*y+y^3*z+z^4;
//f:=y^3*z-(x^4+3*x^3*z+z^4);
//f:=x^3*y+y^3*z+z^3*x;
// f:=y^3*z+z^4-x^4;
f:=x^4+y^4+z^4;
C:=Curve(PP,f);

T:=Twists(C);
#T;
[[[i,j] : j in [i+1..#T], i in [1..#T-1] | IsIsomorphic(T[i],T[j])]];

exit;

// part to be integrated
_,Aut:= IsIsomorphicTernaryQuartics(f,f : geometric:=true);
e:=Lcm([Degree(BaseRing(M)) : M in Aut]);
p:=Characteristic(BaseRing(Aut[1]));
F:=GF(p^e);
for M in Aut do Embed(BaseRing(M),F); end for;
// Be careful about the transpose to have the right action of automorphims
Aut:=[NormalizedM(ChangeRing(Transpose(M),GF(p^e))) : M in Aut];
#Aut;


T:=Twists(C,Aut);
#T;
//[LPolynomial(D) : D in T];
[[[i,j] : j in [i+1..#T], i in [1..#T-1] | IsIsomorphic(T[i],T[j])]];
