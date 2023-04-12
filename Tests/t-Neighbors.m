
//A 2-neighbor step with only An fibers at t=0: Hilbert D=17

A<g,h> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([(1+2*g*t+(2*h+(g+1)^2)*t^2+2*(g*h+g+2*g^2+h)*t^3+((g+h)^2+2*g^3)*t^4),
                 - 4*h^2*t^5*(1 + g*t + (h + 2*g+1)*t^2 + (h+2*g^2+g)*t^3 ),
                 4*h^4*t^10*( (2*g+1)*t^2 + 1)]);
RootConfiguration(S);
a := 2;
v := Frame(S) ! ([-1,-2,-3] cat [-4: i in [1..10]] cat [-3,-2,-1]);

//Elliptic parameter is x/t^4.
u := NewEllipticParameter(S, v, a);
u;
x := Parent(u).1;
assert u eq x/t^4;
