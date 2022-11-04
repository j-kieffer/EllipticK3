
A<r> := PolynomialRing(Rationals());
F<r> := FieldOfFractions(A);
P<t> := PolynomialRing(F);

/* Example with An fibers: Shimura N=6 */

S := EllipticK3(F, t, [(r+1)*t^2 + (3*r^2+2*r)*t + r^3,
		       (r+1)*((4*r+2)*t + 2*r^2) * (t^3-t^2),
		       (r+1)^2*(t+r)*(t^3-t^2)^2]);

S2 := EllipticK3(F, t, [(r+1)^2*t^2 + (3*r^4+4*r^3+2*r^2)*t + r^6,
		       2*(r+1)^2*((2*r^2+2*r+1)*t + r^4)*(t-(2*r+1))*t^2,
		       (r+1)^4*(t-(2*r+1))^2*(t+r^2)*t^4]);

print KodairaType(S2, t);
Pl := t-(2*r+1);
print KodairaType(S2, Pl);

HasRationalComponents(S, t, "I7");
assert HasRationalComponents(S2, t, "I7");
HasRationalComponents(S, t-1, "I4");
assert HasRationalComponents(S2, Pl, "I4");

/* Example with E6 and D5: Hilbert D=28 (typo in c in paper: cf aux files) */

A<f,g> := PolynomialRing(Rationals(), 2);
F<f,g> := FieldOfFractions(A);
P<t> := PolynomialRing(F);
a := 2*(f^2-g^2)*(t-1)+t;
b := (f^2-g^2)^2*(1-t)-2*(f^2-g^2)*(f+1)*t;
c := (f+1)^2*(f^2-g^2)^2*t;
S := EllipticK3(F, t, [a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);

print KodairaType(S, P!0);
print KodairaType(S, t);

assert HasRationalComponents(S, t, "I1*");
S := InvertT(S);
assert HasRationalComponents(S, t, "IV*");

/* Another E6/D5: Hilbert D=37 */

a := (2*g-f+1)*(t-1)+g^2*t/4;
b := (f-g-1)*((f+g-1)*(t-1) + (f-2)*g*t/2);
c := (g-f+1)^2*(f^2*(t-1)+(f-2)^2)/4;
S := EllipticK3(F, t, [a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]);

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
S := EllipticK3(F, t, [a, 2*b*t^4, c*t^8]);

print KodairaType(S, P!0);
print KodairaType(S, t);
assert HasRationalComponents(S, t, "I11");
S := InvertT(S);
assert HasRationalComponents(S, t, "I2*");

