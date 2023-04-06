
/* Example: Hilbert D=77 */

AttachSpec("spec");

//Construct the elliptic K3 surface
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

//Let's add a Mordell--Weil section of height 11/12
X1 := -a3/2*t*(t-1)^2;
AddSection(~S, X1);
MordellWeilRank(S);
Height(S, MordellWeilSections(S)[1]);

//We now add the 2-torsion section
Xtors := PolynomialRing(S)!0;
AddSection(~S, Xtors);
#TorsionGroup(S);
TorsionSections(S);

//And finally the second Mordell--Weil section of height 7/4
X2 := (r*s^2 - r)^2*t^3*(t-1);
AddSection(~S, X2);
MordellWeilRank(S);
Height(S, MordellWeilSections(S)[2]);
NeronSeveriDiscriminant(S);

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

//Another example: the CM point of discriminant -627.

K := QuadraticField(-11); //So that both sections are defined.
AssignNames(~K, ["sqrt11"]);
P<t> := PolynomialRing(K);
r := -7/4;
p := 4*(r-1)*(r^2-2)+1;
d := (r^2-1)^2*(9*t+(2*r-1)*p);
c := 9*t^2-(2*r-1)*(8*r^2+4*r-22)*t + (2*r-1)^2*p;
b := (t-(r^2-2*r))*c+d;
a := (t-(r^2-2*r))^2*c + 2*(t-(r^2-2*r))*d + (r^2-1)^4*((4*r+4)*t+p);
S := EllipticK3([a, 8*(r-1)^4*(r+1)^5*b*t^2, 16*(r-1)^8*(r+1)^10*c*t^4]);

X1 := -4*(r-1)^4*(r+1)^5*(2*r-1)*t^2/(r^2-r+1)^2
      + 4*(r-2)*(r+1)^4*t^3/(r^2-r+1);
Y1 := Sqrt(RHS(S,X1));
AddSection(~S, X1, Y1);

q := 419430400*t^5 + 2846883840*t^4 + 17148174336*t^3 + 78784560576*t^2 +
     175272616341*t - 12882888;
X2 := 3^5*11^4*t^2*q/(2^12*(81920*t^3 + 9216*t^2 + 23868*t + 39339)^2);
Y2 := Sqrt(RHS(S,X2));
AddSection(~S, X2, Y2);

MordellWeilRank(S);
Height(S, MordellWeilSections(S)[1]);
Height(S, MordellWeilSections(S)[2]);
TorsionGroup(S);
NeronSeveriDiscriminant(S);
