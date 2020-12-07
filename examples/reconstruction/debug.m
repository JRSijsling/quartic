FF := FiniteField(23);
S<x1,x2,x3> := PolynomialRing(FF, 3);

f := 10*x1^4 + 7*x1^3*x3 + 13*x1^2*x2^2 + 2*x1^2*x2*x3 + 17*x1^2*x3^2 + 11*x1*x2^3 + 22*x1*x2^2*x3 + 11*x1*x2*x3^2 + 11*x1*x3^3 + 2*x2^4 + 20*x2^3*x3 + 19*x2^2*x3^2 + 5*x2*x3^3 + 9*x3^4;
//I := DixmierOhnoInvariants(f);
PP2 := ProjectiveSpace(S);
X := Curve(PP2, f);
I := DixmierOhnoInvariants(X);

print I;
print PlaneQuarticFromDixmierOhnoInvariants(I);
