

intrinsic aInvariants(S :: EllK3) -> SeqEnum[RngUPolElt]

{Return coefficients of the elliptic K3 surface}
    
    F := S`BaseField;
    a2,a4,a6 := Explode(S`Coeffs);
    return [F!0, a2, F!0, a4, a6];

end intrinsic;


intrinsic EllipticCurve(S :: EllK3) -> CrvEll

{Return the generic fiber of S as an elliptic curve}

    R := Parent(S`EllParam);
    Q := FieldOfFractions(R);
    AssignNames(~Q, [Sprint(S`EllParam)]);
    Coeffs := [Q!c : c in aInvariants(S)];
    return EllipticCurve(Coeffs);

end intrinsic;


intrinsic cInvariants(S :: EllK3) -> SeqEnum[RngUPolElt]
								
{c4, c6 of the elliptic K3 surface}
    
    //Get c4,c6 as formal polynomials
    R<a2,a4,a6> := PolynomialRing(Integers(), 3);
    F := FieldOfFractions(R);
    E := EllipticCurve([F!x: x in [0,a2,0,a4,a6]]);
    FormalC := cInvariants(E);

    return [Evaluate(R!Numerator(c), S`Coeffs): c in FormalC];
    
end intrinsic;


intrinsic Discriminant(S :: EllK3) -> RngUPolElt
							 
{Discriminant of the elliptic K3 surface}
    
    //Get discriminant as a formal polynomial
    R<a2,a4,a6> := PolynomialRing(Integers(), 3);
    F := FieldOfFractions(R);
    E := EllipticCurve([F!x: x in [0,a2,0,a4,a6]]);
    D := R!Numerator(Discriminant(E));

    return Evaluate(D, S`Coeffs);

end intrinsic;
