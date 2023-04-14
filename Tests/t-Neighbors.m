
//A 2-neighbor step with only An fibers at t=0: Hilbert D=17

A<g,h> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([(1+2*g*t+(2*h+(g+1)^2)*t^2+2*(g*h+g+2*g^2+h)*t^3+((g+h)^2+2*g^3)*t^4),
                 - 4*h^2*t^5*(1 + g*t + (h + 2*g+1)*t^2 + (h+2*g^2+g)*t^3 ),
                 4*h^4*t^10*( (2*g+1)*t^2 + 1)]);
v := Frame(S) ! ([-1,-2,-3] cat [-4: i in [1..10]] cat [-3,-2,-1]);

//Elliptic parameter is x/t^4.
u := NewEllipticParameter(S, v);
x := Parent(u).1;
assert u eq x/t^4;
u;

//Hilbert D=8: D9 + E7, identifying E8 fiber

A<r,s> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([t*( (2*r+1)*t+r ), 2*r*s*t^4*(t+1), r*s^2*t^7]);
v := Frame(S) ! ([0,0,0,0,0,0,0,-4,-4,-7,-6,-5,-4,-3,-2,-1]);

u := NewEllipticParameter(S, v);
u;

//Hilbert D=12: E8+D6+A2, identifying E7 across 2 fibers

A<e,f> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([((1-f^2)*(1-t) + t)*t,
                 2*e*t^3*(t-1),
                 e^2*(t-1)^2*t^5]);

v := Frame(S) ! [0,0,0,0,0,0,0,0,-3,-2,-4,-3,-2,-1,-1,-1];
u := NewEllipticParameter(S,v);
u;

//Hilbert D=28: E6+D5+A4+section, identifying D8 that doesn't touch the section

A<f,g> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
a := 2*(f^2-g^2)*(t-1)+t;
b := (f^2-g^2)^2*(1-t)-2*(f^2-g^2)*(f+1)*t;
c := (f+1)^2*(f^2-g^2)^2*t;
S := EllipticK3([a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);

X := P!0;
Y := Sqrt(RHS(S,X));
AddSection(~S, X, Y);
v := Frame(S) ! [0,0,0,0,0,0,-1,-1,-2,-2,-2,-1,-2,-2,-1,0];
u := NewEllipticParameter(S,v);
u;

//Hilbert D=37: same as 28

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
AddSection(~S, X, Y);

v := Frame(S) ! [0,0,0,0,0,0,-1,-1,-2,-2,-2,-1,-2,-2,-1,0];
u := NewEllipticParameter(S,v);
u;

//Hilbert D=53: A8+A6+A1+section, identifying E7

A<g,h> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);

a := g^2*t^4 + (4*(h+1)^2 - 2*g)*t^3 + (4*h^2+4*g*h+6*g-3)*t^2 + 2*(8*h^2-4*g*h+8*h+1)*t + (2*h+1)^2;
b := -4*(h-g+1)*( g*(2*h+1)*t^2 + (6*h^2-2*g*h+6*h+1)*t + (2*h + 1)^2);
c := 16*(h-g+1)^2*(2*h+1)^2;
S := EllipticK3([a, 2*b*t*(t-h), c*t^2*(t-h)^2]);
X := P!0;
Y := Sqrt(RHS(S,X));
AddSection(~S, X, Y);
v := Frame(S) ! [-1,-2,-3,-4,-4,-3,-2,-1,0,0,0,0,0,0,0,0];
u := NewEllipticParameter(S, v);
u;
