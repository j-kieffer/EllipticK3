
P<a,b> := PolynomialRing(Rationals(), 2);
Q := FieldOfFractions(P);
AssignNames(~Q, ["a","b"]);
R<t> := PolynomialRing(Q);

S := EllipticK3(Q, t, [a*t^2, b*t^7, (a+b)*t^10]);
print S;

Q<t> := PolynomialRing(Rationals());

S := EllipticK3(Rationals(), t, [t^5, t^9]);
print S;

