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
