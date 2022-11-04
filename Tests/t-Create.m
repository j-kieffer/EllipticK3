
P<a,b> := PolynomialRing(Rationals(), 2);
R<t> := PolynomialRing(P);

S := EllipticK3([a*t^2, b*t^7, (a+b)*t^10]);
assert IsField(S`BaseField);

Q<t> := PolynomialRing(Rationals());

S := EllipticK3([t^5, t^9]);
assert S`EllParam eq t;
