
//A 2-neighbor step with only An fibers at t=0: Hilbert D=17

A<g,h> := RationalFunctionField(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([(1+2*g*t+(2*h+(g+1)^2)*t^2+2*(g*h+g+2*g^2+h)*t^3+((g+h)^2+2*g^3)*t^4),
                 - 4*h^2*t^5*(1 + g*t + (h + 2*g+1)*t^2 + (h+2*g^2+g)*t^3 ),
                 4*h^4*t^10*( (2*g+1)*t^2 + 1)]);
assert NeronSeveriDiscriminant(S) eq 17;
v := Frame(S) ! ([-1,-2,-3] cat [-4: i in [1..10]] cat [-3,-2,-1]);
u := NewEllipticParameter(S, v);
x := Parent(u).1;
uu := x/t^4;
//The following asserts that u and uu are indeed scalar multiples.
assert Degree(PolynomialRing(S) ! (u/uu)) eq 0;

//Hilbert D=8: D9 + E7, identifying E8 fiber

A<r,s> := RationalFunctionField(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([t*( (2*r+1)*t+r ), 2*r*s*t^4*(t+1), r*s^2*t^7]);
assert NeronSeveriDiscriminant(S) eq 8;
v := Frame(S) ! ([0,0,0,0,0,0,0,-4,-4,-7,-6,-5,-4,-3,-2,-1]);
u := NewEllipticParameter(S, v);
x := Parent(u).1;
uu := (x+s*t^3)/t^4;
assert Degree(PolynomialRing(S) ! (u/uu)) eq 0;

//Hilbert D=12: E8+D6+A2, identifying E7 across 2 fibers

A<e,f> := RationalFunctionField(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([((1-f^2)*(1-t) + t)*t,
                 2*e*t^3*(t-1),
                 e^2*(t-1)^2*t^5]);
assert NeronSeveriDiscriminant(S) eq 12;
v := Frame(S) ! [0,0,0,0,0,0,0,0,-2,-3,-4,-3,-2,-1,-1,-1];
//Elliptic parameter is 
u := NewEllipticParameter(S,v);
x := Parent(u).1;
uu := (x + e/(f+1)*t^2*(t-1))/(t^3*(t-1));
assert Degree(PolynomialRing(S) ! (u/uu)) eq 0;

//Hilbert D=28: E6+D5+A4+section, identifying D8 that doesn't touch the section

A<f,g> := RationalFunctionField(Rationals(), 2);
P<t> := PolynomialRing(A);
a := 2*(f^2-g^2)*(t-1)+t;
b := (f^2-g^2)^2*(1-t)-2*(f^2-g^2)*(f+1)*t;
c := (f+1)^2*(f^2-g^2)^2*t;
S := EllipticK3([a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);
X := P!0;
Y := Sqrt(RHS(S,X));
AddSection(~S, X, Y);
assert NeronSeveriDiscriminant(S) eq 28;
v := Frame(S) ! [0,0,0,0,0,0,-1,-1,-2,-2,-2,-1,-2,-2,-1,0];
u := NewEllipticParameter(S,v);
x := Parent(u).1;
uu := x/(t^2*(t-1)^2);
assert Degree(PolynomialRing(S) ! (u/uu)) eq 0;

//Hilbert D=37: same as 28

A<f,g> := RationalFunctionField(Rationals(), 2);
P<t> := PolynomialRing(A);
a1 := g^2/4;
r0 := g - f + 1;
s0 := f -1;
r1 := (f-2)*r0/g;
a := (t-1)*(2*r0+s0) + a1*t;
b := (1-t)*r0*(r0+2*s0) - 2*a1*r1*t;
c := (t-1)*r0^2*s0 + a1*r1^2*t;
S := EllipticK3([a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);
X := r0*t*(1-t);
Y := Sqrt(RHS(S,X));
AddSection(~S, X, Y);
assert NeronSeveriDiscriminant(S) eq 37;
v := Frame(S) ! [0,0,0,0,0,0,-1,-1,-2,-2,-2,-1,-2,-2,-1,0];
u := NewEllipticParameter(S,v);
x := Parent(u).1;
uu := (x - r0*t*(t-1)^2 + (g-2*f+2)*t*(t-1)^2)/(t^2*(t-1)^2);
assert Degree(PolynomialRing(S) ! (u/uu)) eq 0;

//Hilbert D=53: A8+A6+A1+section, identifying E7

A<g,h> := RationalFunctionField(Rationals(), 2);
P<t> := PolynomialRing(A);
a := g^2*t^4 + (4*(h+1)^2 - 2*g)*t^3 + (4*h^2+4*g*h+6*g-3)*t^2 + 2*(8*h^2-4*g*h+8*h+1)*t + (2*h+1)^2;
b := -4*(h-g+1)*( g*(2*h+1)*t^2 + (6*h^2-2*g*h+6*h+1)*t + (2*h + 1)^2);
c := 16*(h-g+1)^2*(2*h+1)^2;
S := EllipticK3([a, 2*b*t*(t-h), c*t^2*(t-h)^2]);
X := P!0;
Y := Sqrt(RHS(S,X));
AddSection(~S, X, Y);
assert NeronSeveriDiscriminant(S) eq 53;
v := Frame(S) ! [-1,-2,-3,-4,-4,-3,-2,-1,0,0,0,0,0,0,0,0];
u := NewEllipticParameter(S, v);
x := Parent(u).1;
uu := x;
assert Degree(PolynomialRing(S) ! (u/uu)) eq 0;

//Hilbert D=61: D7+A6+A2+section, identifying E7

A<g,h> := RationalFunctionField(Rationals(), 2);
P<t> := PolynomialRing(A);
c := 16*(g+1)^2*(h-g)^4*(g*h*t+2*h*t-g^2*t-2*g*t+g^2+g)^2;
b2 := -4*(g+1)*(h-g)^4*(2*g*h^2+4*h^2+g^2*h-g^2-6*g-6);
b1 := -4*(g+1)^2*(h-g)^3*(g^2*h+2*h-2*g^2-6*g);
b0 := 4*g^2*(g+1)^3*(h-g)^2;
a3 := 4*h^3*(h-g)^3;
a2 := (h-g)^2*(g^2*h^2-4*g*h^2-8*h^2-2*g^2*h+g^2+12*g+12);
a1 :=  -2*(g+1)*(h-g)*(g^2*h+4*h-g^2-6*g);
a0 := g^2*(g+1)^2;
S := EllipticK3([(a0 + a1*t + a2*t^2 + a3*t^3),
                 2*(b0 + b1*t + b2*t^2)*t*(t-1),
                 c*t^2*(t-1)^2]);
X := P!0;
Y := Sqrt(RHS(S,X));
AddSection(~S, X, Y);
assert NeronSeveriDiscriminant(S) eq 61;
v := Frame(S) ! [-3,-3,-5,-4,-3,-2,-1,0,0,0,0,0,0,-1,-1,0];
u := NewEllipticParameter(S, v);
x := Parent(u).1;
uu := x/(t-1);
assert Degree(PolynomialRing(S) ! (u/uu)) eq 0;

//Hilbert D=65: E7+A4+A4+section, identifying E8

A<r,s> := RationalFunctionField(Rationals(), 2);
P<t> := PolynomialRing(A);
v := (s^4 + 2*s^2 - 3)*r + (s^3 + s)/2;
g := (2*s*v-s^2-1)/(s^2-1);
f := (s^2-5)*(s^2*v+v-2*s)/((s-1)*(s+1)*(s^2+3));
b1 := (v^2 - f^2 + g^2 + 1)/2;
e := (f-v-g-1)*(f-v+g+1);
a1 := 2*b1 - 1 - e - (v-f)^2;
c := e^2*(v^2*(1-t) + t);
b := e*(v*f*(1-t) + b1*t*(1-t) + t^2);
a := f^2*(1-t) + a1*t*(1-t) + t^2;
S := EllipticK3([a,2*b*t^2*(t-1),c*t^4*(t-1)^2]);
k := (s-1)*(s+1)*(2*r*s^2-4*r*s+s+2*r-2)*(2*r*s^2+4*r*s+s+2*r+2)
     *(2*r*s^3-2*r*s^2+s^2+6*r*s-s-6*r+2)
     *(2*r*s^3+2*r*s^2+s^2+6*r*s+s+6*r+2)
     /(16*(2*r*s^2+s+2*r)^2);
X := k*t^2;
Y := Sqrt(RHS(S,X));
AddSection(~S, X, Y);
assert NeronSeveriDiscriminant(S) eq 65;
v := Frame(S) ! [-3,-4,-5,-6,-4,-2,-3,0,0,0,0,-1,-1,-1,-1,0];
u := NewEllipticParameter(S, v);
x := Parent(u).1;
//Since we choose to kill the t coefficient, elliptic parameter is not exactly this:
//uu := (x-k*t*(t-1))/(t-1);
uu := (x - k*t*(t-1))/(t-1) - k;
assert Degree(PolynomialRing(S) ! (u/uu)) eq 0;

//This is the other E8.
//v := Frame(S) ! [-3,-4,-5,-6,-4,-2,-3,-1,-1,-1,-1,0,0,0,0,0];
//u := NewEllipticParameter(S, v);

//Now look at examples with sections

//Hilbert D=13: E8+E6+A1+section, identifying E7 on section + 1 fiber

A<g,h> := RationalFunctionField(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([(4*g+1)*t^2, -4*g*(h-g-1)*(t-1)*t^3, 4*g^2*(t-1)^2*t^4*(h^2*t+1)]);
X := (t-1)*t^2*((2-h)*h + h^2*t);
Y := t^2*(t-1)*(h^3*t^3 -h^2*(2*h-3)*t^2 + h*(h^2-3*h+2*g+2)*t + 2*g);
AddSection(~S, X, Y);
v := Frame(S) ! [0,0,0,0,0,0,0,0,-2,-2,-2,-1,0,-1,0,1];
