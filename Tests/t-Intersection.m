
/* Hilbert D=5 */

A<g,h> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
a := P!0;
b := 1/4*t^3*(-3*g^2*t+4);
c := -1/4*t^5*(4*h^2*t^2 + (4*h+g^3)*t + (4*g+1));
S := EllipticK3([a, b, c]);

assert RootConfiguration(S) eq RootLatticeType("E7") + RootLatticeType("E8");
X := 1/4*t^2*((1+2*h*t)^2 + 4*g);
Y := 1/8*t^3*(1+2*h*t)*((1+2*h*t)^2+6*g);
assert Y^2 eq RHS(S,X);

F0 := ReducibleFiber(S, t);
assert Intersection(S, F0, X, Y) eq 1;
assert IsPositiveDefinite(IntersectionMatrix(S,X,Y));

/* Hilbert D=28 */

A<f,g> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
a := 2*(f^2-g^2)*(t-1)+t;
b := (f^2-g^2)^2*(1-t)-2*(f^2-g^2)*(f+1)*t;
c := (f+1)^2*(f^2-g^2)^2*t;
S := EllipticK3([a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);

X := P!0;
Y := Sqrt(RHS(S,X));
F, F0, F1 := Explode(ReducibleFibers(S));
assert Intersection(S, F, X, Y) in [1,2];
assert Intersection(S, F0, X, Y) eq 1;
assert Intersection(S, F1, X, Y) in [2,3];

/* Hilbert D=41 */

A<r,s> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);

a := t^4 + ((-2*s+1/2)*r + (16*s+2)) * t^3 + ((s^2+1/2*s+1/16)*r^2 +
  (-24*s^2-4*s+1/2)*r + (96*s^2+16*s+1)) * t^2 + ((8*s^3+4*s^2+1/2*s)*r^2 +
  (-96*s^3-32*s^2+2*s)*r + (256*s^3+32*s^2)) * t + ((16*s^4+8*s^3+s^2)*r^2 +
  (-128*s^4-32*s^3)*r + 256*s^4);
b := 4*s*r*t^3 + ((-8*s^2+2*s)*r^2 + (48*s^2+4*s)*r) * t^2 +
  ((4*s^3+2*s^2+1/4*s)*r^3 + (-64*s^3-12*s^2+s)*r^2 + (192*s^3+16*s^2)*r) * t +
  ((16*s^4+8*s^3+s^2)*r^3 + (-128*s^4-32*s^3)*r^2 + 256*s^4*r);
c := 16*s^2*r^2*t^2 + ((-32*s^3+8*s^2)*r^3 + 128*s^3*r^2) * t +
  ((16*s^4+8*s^3+s^2)*r^4 + (-128*s^4-32*s^3)*r^3 + 256*s^4*r^2);
S := EllipticK3([a, 2*b*t^2, c*t^4]);
X := -4*r*s*t;
Y := Sqrt(RHS(S,X));

F, F0 := Explode(ReducibleFibers(S));
assert Intersection(S,F0,X,Y) in [1,5];
assert Intersection(S,F,X,Y) in [4,7];
