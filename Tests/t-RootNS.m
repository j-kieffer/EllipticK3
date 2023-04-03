
/* Example with E6 and D5: Hilbert D=28 (typo in c in paper: cf aux files) */

A<f,g> := PolynomialRing(Rationals(), 2);
F<f,g> := FieldOfFractions(A);
P<t> := PolynomialRing(F);
a := 2*(f^2-g^2)*(t-1)+t;
b := (f^2-g^2)^2*(1-t)-2*(f^2-g^2)*(f+1)*t;
c := (f+1)^2*(f^2-g^2)^2*t;
S := EllipticK3([a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);

assert KodairaType(S, P!0) eq "IV*";
assert KodairaType(S, t) eq "I1*";
assert HasRationalComponents(S, t, "I1*");
S := InvertT(S);
assert HasRationalComponents(S, t, "IV*");
