
//A 2-neighbor step with only An fibers at t=0: Hilbert D=17

A<g,h> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([(1+2*g*t+(2*h+(g+1)^2)*t^2+2*(g*h+g+2*g^2+h)*t^3+((g+h)^2+2*g^3)*t^4),
                 - 4*h^2*t^5*(1 + g*t + (h + 2*g+1)*t^2 + (h+2*g^2+g)*t^3 ),
                 4*h^4*t^10*( (2*g+1)*t^2 + 1)]);
RootConfiguration(S);
v := Frame(S) ! ([-1,-2,-3] cat [-4: i in [1..10]] cat [-3,-2,-1]);

//Elliptic parameter is x/t^4.
u := NewEllipticParameter(S, v);
u;
x := Parent(u).1;
assert u eq x/t^4;

//Hilbert D=8: D9 + E7, identifying E8 fiber

A<r,s> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([t*( (2*r+1)*t+r ), 2*r*s*t^4*(t+1), r*s^2*t^7]);
RootConfiguration(S);
v := Frame(S) ! ([0,0,0,0,0,0,0,-4,-4,-7,-6,-5,-4,-3,-2,-1]);

u := NewEllipticParameter(S, v);
u;

//Hilbert D=12

A<e,f> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
S := EllipticK3([((1-f^2)*(1-t) + t)*t,
                 2*e*t^3*(t-1),
                 e^2*(t-1)^2*t^5]);
RootConfiguration(S);

v := Frame(S) ! [0,0,0,0,0,0,0,0,-3,-2,-4,-3,-2,-1,-1,-1];
u := NewEllipticParameter(S,v);
u;
