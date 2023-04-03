AttachSpec("spec");
P<b> := PolynomialRing(Rationals());
R<t> := PolynomialRing(P);
S := EllipticK3([t, 2*b*t^3*(t-1), b^2*t^5*(t-1)^2]:
                ComputeReducibleFibers := false);

S;
EllipticCurve(S);
BaseField(S);
BaseFunctionField(S);
EllipticParameter(S);
PolynomialRing(S);
aInvariants(S);
Coefficients(S);
cInvariants(S);
Discriminant(S);
RHS(S, 0); //Evaluate RHS of equation of S at x=0
ReduciblePlaces(S);
