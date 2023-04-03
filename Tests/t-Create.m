
P<a,b> := PolynomialRing(Rationals(), 2);
R<t> := PolynomialRing(P);

S := EllipticK3([a*t^2+1, b*t^7, (a+b)*t^10]:
                ComputeReducibleFibers := false);
assert IsField(BaseField(S));

Q<t> := PolynomialRing(Rationals());

S := EllipticK3([t^5+1, t^9]:
                ComputeReducibleFibers := false);
assert EllipticParameter(S) eq t;
