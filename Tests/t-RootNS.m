
/* Example with E6 and D5: Hilbert D=28 (typo in c in paper: cf aux files) */

A<f,g> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
a := 2*(f^2-g^2)*(t-1)+t;
b := (f^2-g^2)^2*(1-t)-2*(f^2-g^2)*(f+1)*t;
c := (f+1)^2*(f^2-g^2)^2*t;
S := EllipticK3([a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);

assert KodairaType(S, P!0) eq "IV*";
assert KodairaType(S, t) eq "I1*";
assert RootConfiguration(S) eq RootLatticeType("E6") + RootLatticeType("D5") + RootLatticeType("A4");
assert Rank(RootLattice(S)) eq 15;
assert NeronSeveriRank(S) eq 17;
assert #ReducibleFibers(S) eq 3;
assert #TorsionGroup(S) eq 1;
assert MordellWeilRank(S) eq 0;
