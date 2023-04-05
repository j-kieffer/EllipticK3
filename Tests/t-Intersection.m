
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

/* Hilbert D=13 */

A<g,h> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([(4*g+1)*t^2, -4*g*(h-g-1)*(t-1)*t^3, 4*g^2*(t-1)^2*t^4*(h^2*t+1)]);
X := (t-1)*t^2*((2-h)*h + h^2*t);
Y := t^2*(t-1)*(h^3*t^3 -h^2*(2*h-3)*t^2 + h*(h^2-3*h+2*g+2)*t + 2*g);
F, F0, F1 := Explode(ReducibleFibers(S));
assert Intersection(S, F, X, Y) eq 0;
assert Intersection(S, F0, X, Y) in [1,2];
assert Intersection(S, F1, X, Y) eq 1;

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

/* Hilbert D=29 */

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
assert Intersection(S, F, X, Y) eq 1;
assert Intersection(S, F0, X, Y) in [1,8];

/* Hilbert D=37 */

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
assert Intersection(S, F, X, Y) in [1,2];
assert Intersection(S, F0, X, Y) in [2,3];
assert Intersection(S, F1, X, Y) in [1,4];

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

/* Hilbert D=56 */

A<g,h> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);

la := (2*h - g^2 + 1);
mu := -2*h;
c := -(g^2-1)^2*( 4*(g^2-1)*t-(h+2)^2);
b := (g^2-1)*( 2*(g^2-1)*t^2 -(h^2+4*h-4*g^2+8)*t - (h+2)^2 );
a := (g^2-1)^2*t^3+ (h+2)*(h-2*g^2+4)*t^2 + 2*(h^2+4*h-2*g^2+6)*t + (h+2)^2;
S := EllipticK3([a, 2*b*t^2*(la*t-mu), c*t^4*(la*t-mu)^2]);
w := (g^2 - 1)/(2*h) - 1;
s := 2/w;
v := -(g^2-1)*s/((g^2+3)*s+8);
t := EllipticParameter(S);
X := v^2*s*(s+2)*t^2 -v^2*s*(s+2)^2/2*t^3; //Rescaling needed.
//Y := Sqrt(RHS(S,X)); 

/* Hilbert D=77 */

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
assert Intersection(S,F,X1,Y1) in [2,3];
assert Intersection(S,F0,X1,Y1) in [1,5];
assert Intersection(S,F1,X1,Y1) eq 2;
assert Intersection(S,Fa,X1,Y1) eq 0;

/* Hilbert D=88 */

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
assert Intersection(S,F,X,Y) ne 0;
assert Intersection(S,F0,X,Y) in [2,8];
assert Intersection(S,F1,X,Y) in [2,3];
