
A<g,h> := PolynomialRing(Rationals(), 2);
F<g,h> := FieldOfFractions(A);
R<t> := PolynomialRing(F);

S := EllipticK3([1/4*t^3*(-3*g^2*t+4),
		 -1/4*t^5*(4*h^2*t^2 + (4*h+g^3)*t + (4*g+1))]);

x := 1/4*t^2*((1+2*h*t)^2+4*g);
y := 1/8*t^3*(1+2*h*t)*((1+2*h*t)^2+6*g);

assert RHS(S, x) eq y^2;
