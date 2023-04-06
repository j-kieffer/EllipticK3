
/* A section of naive height 6: CM -627 on Shimura N=57, also with non-rat'l fibers */

P<t> := PolynomialRing(Rationals());
r := -7/4;
p := 4*(r-1)*(r^2-2)+1;
d := (r^2-1)^2*(9*t+(2*r-1)*p);
c := 9*t^2-(2*r-1)*(8*r^2+4*r-22)*t + (2*r-1)^2*p;
b := (t-(r^2-2*r))*c+d;
a := (t-(r^2-2*r))^2*c + 2*(t-(r^2-2*r))*d + (r^2-1)^4*((4*r+4)*t+p);
S := EllipticK3([a, 8*(r-1)^4*(r+1)^5*b*t^2, 16*(r-1)^8*(r+1)^10*c*t^4]);

ReducibleFibers(S);
F, F0 := Explode(ReducibleFibers(S));

FieldOfRationality(F);
FieldOfRationality(F0);

Components(F);
Components(F0);

X1 := -4*(r-1)^4*(r+1)^5*(2*r-1)*t^2/(r^2-r+1)^2
      + 4*(r-2)*(r+1)^4*t^3/(r^2-r+1);
Y1 := Sqrt(RHS(S,X1));
assert Intersection(F,X1,Y1) in [1,11];
assert Intersection(F0,X1,Y1) eq 3;
