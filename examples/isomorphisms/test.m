P<x,y,z> := PolynomialRing(Rationals(), 3);
PP := ProjectiveSpace(P);
f := x^4 + 3*y^4 + 5*z^4 + x^2*y*z + x*y*z^2 + x^2*y^2;
C := Curve(PP, f);
IsIsomorphicQuartic(C, C);
