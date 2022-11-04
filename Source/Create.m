

intrinsic EllipticK3(Coeffs :: SeqEnum[RngUPolElt]) -> EllK3
{Define an elliptic K3 surface by the data of its coefficients.}
    

    require #Coeffs eq 2 or #Coeffs eq 3: "Must specify either 2 or 3
    coefficients";

    S := New(EllK3);

    F := BaseRing(Coeffs[1]);
    P := Parent(Coeffs[1]);

    /* Convert base ring to a field if not already the case */

    if not IsField(F) then	
	FF := FieldOfFractions(F);
	t := P.1;
	P := PolynomialRing(FF);
	try
	    names := [Sprint(F.i): i in [1..NumberOfGenerators(F)]];
	    AssignNames(~FF, names);
	    catch e;
	end try;	
	AssignNames(~P, [Sprint(t)]);
	F := FF;
    end if;
    
    S`BaseField := F;
    S`EllParam := P.1;
    S`Coeffs := [P!c: c in Coeffs];
    
    if #Coeffs eq 2 then S`Coeffs := [P!0] cat S`Coeffs; end if;
    
    return S;

end intrinsic;
	

intrinsic EllipticK3(L :: Lat) -> EllK3
{Construct the family of elliptic K3 surfaces polarized by the given
Néron--Severi lattice.}

    error "Not implemented";

end intrinsic;
