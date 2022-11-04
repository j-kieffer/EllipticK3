
P<b> := PolynomialRing(Rationals());
Q<b> := FieldOfFractions(P);
R<t> := PolynomialRing(Q);

S := EllipticK3(Q, t, [t, 2*b*t^3*(t-1), b^2*t^5*(t-1)^2]);

Pl := ReduciblePlaces(S);

assert #Pl eq 3;
assert t in Pl;
assert t-1 in Pl;
assert 0 in Pl;

assert RootLatticeType(S, t) eq RootConfiguration("D7");
assert RootLatticeType(S, R!0) eq RootConfiguration("E8");
assert RootLatticeType(S, t-1) eq RootConfiguration("A2");

