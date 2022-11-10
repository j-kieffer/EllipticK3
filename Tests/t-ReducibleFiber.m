
A<r> := PolynomialRing(Rationals());
F<r> := FieldOfFractions(A);
P<t> := PolynomialRing(F);

/* Test components: Shimura N=6 */

S := EllipticK3([t, 2*r^2*t^3*(t-1), r^4*t^5*(t-1)^2]);
F := ReducibleFiber(S, t);
assert F`RootType eq RootLatticeType("D7");

x,uk,y1 := Explode(F`RatComps[2]);
assert x eq r^2*t^2;
assert uk eq t^3;

x,uk,y2 := Explode(F`RatComps[3]);
assert x eq r^2*t^2;
assert uk eq t^3;

assert y2 eq -y1;
assert y1^2 eq r^6*t^6;

/* Example with An fibers: Shimura N=14 */

S := EllipticK3([(r+1)*t^2 + (3*r^2+2*r)*t + r^3,
		 (r+1)*((4*r+2)*t + 2*r^2) * (t^3-t^2),
		 (r+1)^2*(t+r)*(t^3-t^2)^2]);

S2 := EllipticK3([(r+1)^2*t^2 + (3*r^4+4*r^3+2*r^2)*t + r^6,
		  2*(r+1)^2*((2*r^2+2*r+1)*t + r^4)*(t-(2*r+1))*t^2,
		  (r+1)^4*(t-(2*r+1))^2*(t+r^2)*t^4]);

assert KodairaType(S2, t) eq "I7";
Pl := t-(2*r+1);
assert KodairaType(S2, Pl) eq "I4";

assert not HasRationalComponents(S, t, "I7");
assert HasRationalComponents(S2, t, "I7");
//print ReducibleFiber(S2,t)`RatComps;
assert not HasRationalComponents(S, t-1, "I4");
assert HasRationalComponents(S2, Pl, "I4");
//print ReducibleFiber(S2,Pl)`RatComps;

/* Example with E6 and D5: Hilbert D=28 (typo in c in paper: cf aux files) */

A<f,g> := PolynomialRing(Rationals(), 2);
F<f,g> := FieldOfFractions(A);
P<t> := PolynomialRing(F);
a := 2*(f^2-g^2)*(t-1)+t;
b := (f^2-g^2)^2*(1-t)-2*(f^2-g^2)*(f+1)*t;
c := (f+1)^2*(f^2-g^2)^2*t;
S := EllipticK3([a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);

assert KodairaType(S, P!0) eq "IV*";
assert KodairaType(S, t) eq "I1*";

assert HasRationalComponents(S, t, "I1*");
S := InvertT(S);
assert HasRationalComponents(S, t, "IV*");

/* Another E6/D5: Hilbert D=37 */

a := (2*g-f+1)*(t-1)+g^2*t/4;
b := (f-g-1)*((f+g-1)*(t-1) + (f-2)*g*t/2);
c := (g-f+1)^2*(f^2*(t-1)+(f-2)^2)/4;
S := EllipticK3([a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);

assert HasRationalComponents(S, t, "I1*");
S := InvertT(S);
assert HasRationalComponents(S, t, "IV*");

/* Dn for n even: Hilbert D=44 */

A<r,s> := PolynomialRing(Rationals(), 2);
F<r,s> := FieldOfFractions(A);
P<t> := PolynomialRing(F);
a := s*(r^2*s^2-s^2+2*s+2)*t^3 - s*(2*r^2*s-3*s+2)*t^2 + s*(r^2-3)*t + 1;
b := s^2*((2*r^2*s^2-2*s^2+4*s+1)/2*t^2 - (r^2*s-2*s+1)*t - 1);
c := s^4*((r^2*s-s+2)*t + 1);
S := EllipticK3([a, 2*b*t^4, c*t^8]);

assert KodairaType(S, P!0) eq "I2*";
assert KodairaType(S, t) eq "I11";
assert HasRationalComponents(S, t, "I11");
S := InvertT(S);
assert HasRationalComponents(S, t, "I2*");

