
P<a,b> := PolynomialRing(Rationals(), 2);
I := Ideal(a^2-b);
Q := quo<P|I>;
F<a,b> := FieldOfFractions(Q);
R<t> := PolynomialRing(F);

S := EllipticK3([1+a*t^2, b*t^7, t^10]: ComputeReducibleFibers:=false);
E := EllipticCurve(S);

assert aInvariants(S) eq aInvariants(E);
assert cInvariants(S) eq cInvariants(E);
assert Discriminant(S) eq Discriminant(E);
