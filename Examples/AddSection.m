
/* Example: Hilbert D=77 */

AttachSpec("spec");

//Construct the elliptic K3 surfaces
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

//At this point, S only knows about its root lattice
RootConfiguration(S);
Determinant(RootLattice(S));

//Let's add a Mordell--Weil section:
X1 := -a3/2*t*(t-1)^2;
Y1 := Sqrt(RHS(S,X1));
AddSection(~S, X1, Y1);

//We now add the other 2-torsion section
Xtors := PolynomialRing(S)!0;
Ytors := Sqrt(RHS(S,Xtors));
AddSection(~S, Xtors, Ytors);

//We can get trivial info about the section:
xy := MordellWeilSections(S)[1];
xy;
NSVector(xy);
xCoordinate(xy);
yCoordinate(xy);
GenericPoint(xy);

//And also some less trivial info
Height(S, xy);
MordellWeilRank(S);
NeronSeveriRank(S);
TorsionGroup(S);
MordellWeilGroup(S);
Frame(S);
MordellWeilLattice(S);

//We now add the other 2-torsion section
Xtors := PolynomialRing(S)!0;
Ytors := Sqrt(RHS(S,Xtors));
AddSection(~S, Xtors, Ytors);
