

intrinsic EllipticK3(Base :: Rng,
		     EllParam :: RngUPolElt,
		     Coeffs :: SeqEnum[RngUPolElt]) -> EllK3
{Define an elliptic K3 surface by the data of its coefficients.}
    
    require IsField(Base): "Base ring must be a field";

    require #Coeffs eq 2 or #Coeffs eq 3: "Must specify either 2 or 3
    coefficients";

    require BaseRing(EllParam) eq Base and Parent(EllParam).1 eq EllParam:
    "Elliptic parameter must be a variable over the base field";

    S := New(EllK3);
    S`BaseField := Base;
    S`EllParam := EllParam;
    P := Parent(EllParam);

    S`Coeffs := [P!c: c in Coeffs];
    if #Coeffs eq 2 then S`Coeffs := [P!0] cat S`Coeffs; end if;
    
    return S;

end intrinsic;
	

intrinsic EllipticK3(L :: Lat) -> EllK3
{Construct the family of elliptic K3 surfaces polarized by the given
Néron--Severi lattice.}

    error "Not implemented";

end intrinsic;
