
/* Example with E6 and D5: Hilbert D=28 (typo in c in paper: cf aux files) */

AttachSpec("spec");
A<f,g> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
a := 2*(f^2-g^2)*(t-1)+t;
b := (f^2-g^2)^2*(1-t)-2*(f^2-g^2)*(f+1)*t;
c := (f+1)^2*(f^2-g^2)^2*t;
S := EllipticK3([a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);

ReducibleFibers(S);
RootConfiguration(S);
RootLattice(S);
NeronSeveriRootSpan(S) eq RootLattice(S);
Frame(S) eq RootLattice(S);
FrameSpace(S);
NeronSeveriRank(S);
#TorsionGroup(S);
#MordellWeilGroup(S);
MordellWeilRank(S);
MordellWeilSections(S);
