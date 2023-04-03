/* Example with E6 and D5: Hilbert D=28 (typo in c in paper: cf aux files) */

AttachSpec("spec");
A<f,g> := PolynomialRing(Rationals(), 2);
P<t> := PolynomialRing(A);
a := 2*(f^2-g^2)*(t-1)+t;
b := (f^2-g^2)^2*(1-t)-2*(f^2-g^2)*(f+1)*t;
c := (f+1)^2*(f^2-g^2)^2*t;
S := EllipticK3([a*t, b*t^2*(t-1)^2, c*t^3*(t-1)^4]:
                ComputeReducibleFibers := false);
ReduciblePlaces(S);

F := ReducibleFiber(S, P!0);
F;
Place(F);
KodairaType(F);
FieldOfRationality(F);
RootType(F);
RootLattice(F);
Components(F);
Component(F,1);

F0 := ReducibleFiber(S,t);
RootType(F0);
F1 := ReducibleFiber(S,t-1);
RootType(F1);
