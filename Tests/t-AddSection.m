
/* Hilbert D=5 */

A<g,h> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
a := P!0;
b := 1/4*t^3*(-3*g^2*t+4);
c := -1/4*t^5*(4*h^2*t^2 + (4*h+g^3)*t + (4*g+1));
S := EllipticK3([a, b, c]);
X := 1/4*t^2*((1+2*h*t)^2 + 4*g);
Y := 1/8*t^3*(1+2*h*t)*((1+2*h*t)^2+6*g);

M := IntersectionMatrix(S, X, Y);
assert IsPositiveDefinite(M);
ExtendFrameSpace(~S, M);

B := Basis(RootLattice(S));
for v in B do
    assert (FrameSpace(S)!v) * MordellWeilProjection(S) eq 0;
end for;
I := Intersections(S, X, Y);
v := NSVector(S, I);
assert v eq FrameSpace(S) ! ([0: i in [1..15]] cat [1]);

new_frame := ext < Frame(S) | [v] >;
new_rootspan := DualBasisLattice(RootLattice(S)) meet new_frame;
assert new_rootspan eq NeronSeveriRootSpan(S);

AddSection(~S, v, X, Y);
assert NeronSeveriRank(S) eq 18;
assert #TorsionGroup(S) eq 1;
assert MordellWeilRank(S) eq 1;
assert Rank(MordellWeilLattice(S)) eq 1;
assert #MordellWeilSections(S) eq 1;

xy := MordellWeilSections(S)[1];
assert NSVector(xy) in Frame(S);
assert Section(xy) in EllipticCurve(S);
assert xCoordinate(xy) eq BaseFunctionField(S) ! X;
assert yCoordinate(xy) eq BaseFunctionField(S) ! Y;
assert Height(S, xy) eq 5/2;
assert NeronSeveriDiscriminant(S) eq 5;
assert Determinant(NeronSeveriRootSpan(S)) eq 2;

/* Two more complicated sections: Hilbert D=77 */

A<r,s> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);

b0 := (r*s-1)^2 - (s-5*r)^2;
b1 := 16*r^2;
b := r^2*(s^2-1)^2*t^3*(t-1)^2*(b0 + b1*t);
a0 := (r*s-1)^2;
a1 := -(s^2-1)*(r*s - 4*r^2-1);
a2 := r*(s^2-1)*(r*s^2 + 8*s - 57*r)/4;
a3 := 8*r^2*(s^2-1);
a := a0 + a1*t + a2*t^2 + a3*t^3;
S := EllipticK3([a,b,P!0]);
X1 := -a3/2*t*(t-1)^2;
Y1 := Sqrt(RHS(S,X1));
AddSection(~S, X1, Y1);

assert MordellWeilRank(S) eq 1;
assert #MordellWeilSections(S) eq 1;
xy1 := MordellWeilSections(S)[1];
assert Height(S, xy1) eq 11/12;

//Add second section
