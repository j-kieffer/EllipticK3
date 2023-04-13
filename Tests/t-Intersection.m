
/* Hilbert D=5: E7 + E8 + section */

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
assert Intersection(F0, X, Y) eq 1;
assert IsPositiveDefinite(IntersectionMatrix(S,X,Y));

/* Hilbert D=13: E8 + E6 + A1 + section */

A<g,h> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([(4*g+1)*t^2, -4*g*(h-g-1)*(t-1)*t^3, 4*g^2*(t-1)^2*t^4*(h^2*t+1)]);
X := (t-1)*t^2*((2-h)*h + h^2*t);
Y := t^2*(t-1)*(h^3*t^3 -h^2*(2*h-3)*t^2 + h*(h^2-3*h+2*g+2)*t + 2*g);
F, F0, F1 := Explode(ReducibleFibers(S));
assert Intersection(F, X, Y) eq 0;
assert Intersection(F0, X, Y) in [1,5];
assert Intersection(F1, X, Y) eq 1;

/* Hilbert D=28: E6 + D5 + A4 + section */

A<f,g> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
a := 2*(f^2-g^2)*(t-1)+t;
b := (f^2-g^2)^2*(1-t)-2*(f^2-g^2)*(f+1)*t;
c := (f+1)^2*(f^2-g^2)^2*t;
S := EllipticK3([a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);

X := P!0;
Y := Sqrt(RHS(S,X));
F, F0, F1 := Explode(ReducibleFibers(S));
assert Intersection(F, X, Y) in [1,5];
assert Intersection(F0, X, Y) eq 5;
assert Intersection(F1, X, Y) in [2,3];

/* Hilbert D=29: E7 + A8 + section */

A<f,g> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
a := (-(4*f-1)*t^2+(g-2)*t+1);
b := -2*g*t^3*(2*f^2*t^2+(-g+2*f+1)*t-1);
c := g^2*t^6*( (g-4*f)*t+1);
S := EllipticK3([a,b,c]);
t := EllipticParameter(S);
X := t*(-(f-1)*g/f + g*(g-4*f)/(4*f^2)*t);
Y := Sqrt(RHS(S,X));
F, F0 := Explode(ReducibleFibers(S));
assert Intersection(F, X, Y) eq 1;
assert Intersection(F0, X, Y) in [1,8];

/* Hilbert D=37: E6 + D5 + A4 + section */

A<f,g> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
a1 := g^2/4;
r0 := g - f*g - 1;
s0 := f*g + 1;
r1 := f*r0;
a := (t-1)*(2*r0+s0) + a1*t;
b := (1-t)*r0*(r0+2*s0) - 2*a1*r1*t;
c := (t-1)*r0^2*s0 + a1*r1^2*t;
S := EllipticK3([a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);
X := r0*t*(1-t);
Y := Sqrt(RHS(S,X));
F, F0, F1 := Explode(ReducibleFibers(S));
assert Intersection(F, X, Y) in [1,5];
assert Intersection(F0, X, Y) in [1,2];
assert Intersection(F1, X, Y) in [1,4];

/* Hilbert D=41: A5 + A10 + section */

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
assert Intersection(F0,X,Y) in [1,5];
assert Intersection(F,X,Y) in [4,7];

/* Hilbert D=77: A1 + A3 + A5 + D5 + 2-torsion + 2 sections */

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
// Compute second section X2, Y2?
F, F0, F1, Fa := Explode(ReducibleFibers(S));
assert Intersection(F,X1,Y1) in [1,2];
assert Intersection(F0,X1,Y1) in [1,5];
assert Intersection(F1,X1,Y1) eq 2;
assert Intersection(Fa,X1,Y1) eq 0;

/* Hilbert D=88: A9 + D4 + A2 + section */

A<r,s> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);

la := (r+2)^2*s*(2*s+1)^2;
mu := r*s + 2*s + 1;
c1 := 8*r*s^2+16*s^2+8*r*s+8*s+r^2;
c0 := -2;
b3 := 2*r*(r+2)*s*(2*s+1)*(8*s+r^2)*(8*s^2+r);
b2 := 4*s*(8*r^2*s^3+16*r*s^3+24*r^2*s^2+32*r*s^2+16*s^2-r^3*s+20*r^2*s+24*r*s+8*s+r^3-r^2);
b1 := -(16*r*s^2+16*s^2+24*r*s+16*s+r^2);
b0 := 2;
a3 := 4*r*s *(64*r*s^4+16*r^3*s^3+48*r^2*s^3+192*r*s^3+64*s^3+12*r^3*s^2
              +48*r*s^2+r^4*s+12*r^3*s+12*r^2*s+16*r*s+r^3);
a2 := 4*s*(4*r^2*s^3+16*r^2*s^2-8*r*s^2-4*r^3*s+16*r^2*s+16*r*s+4*s-r^3-4*r^2);
a1 := -8*s*(r*s+2*r+1);
a0 := 1;
S := EllipticK3([a0+a1*t+a2*t^2+a3*t^3,
                 8*r^2*s*(b0+b1*t+b2*t^2+b3*t^3)*t^2*(la*t-mu),
                 16*r^4*s^2*t^4*(la*t-mu)^2*(c0+c1*t)^2]);
F, F0, F1 := Explode(ReducibleFibers(S));
X := P!0;
Y := Sqrt(RHS(S,X));
assert Intersection(F,X,Y) in [1,2,4];
assert Intersection(F0,X,Y) in [2,8];
assert Intersection(F1,X,Y) in [2,3];

/* A section of naive height 6: CM -657 on Shimura N=57, also with non-rat'l fibers 
   This is A5 + A11 + 2 sections */

P<t> := PolynomialRing(Rationals());
r := -7/4;
p := 4*(r-1)*(r^2-2)+1;
d := (r^2-1)^2*(9*t+(2*r-1)*p);
c := 9*t^2-(2*r-1)*(8*r^2+4*r-22)*t + (2*r-1)^2*p;
b := (t-(r^2-2*r))*c+d;
a := (t-(r^2-2*r))^2*c + 2*(t-(r^2-2*r))*d + (r^2-1)^4*((4*r+4)*t+p);
S := EllipticK3([a, 8*(r-1)^4*(r+1)^5*b*t^2, 16*(r-1)^8*(r+1)^10*c*t^4]);
F, F0 := Explode(ReducibleFibers(S));
X1 := -4*(r-1)^4*(r+1)^5*(2*r-1)*t^2/(r^2-r+1)^2
      + 4*(r-2)*(r+1)^4*t^3/(r^2-r+1);
Y1 := Sqrt(RHS(S,X1));
assert Intersection(F,X1,Y1) in [1,11];
assert Intersection(F0,X1,Y1) eq 3;

q := 419430400*t^5 + 2846883840*t^4 + 17148174336*t^3 + 78784560576*t^2 +
     175272616341*t - 12882888;
X2 := 3^5*11^4*t^2*q/(2^12*(81920*t^3 + 9216*t^2 + 23868*t + 39339)^2);
Y2 := Sqrt(-11*RHS(S,X2)); //Have to work over Q(sqrt(-11))
assert Intersection(F, X2, Y2) eq 6;
assert Intersection(F0, X2, Y2) eq 3;
