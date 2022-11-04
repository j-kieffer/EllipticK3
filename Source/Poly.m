

intrinsic DegreeOneRoot(pol :: RngUPolElt) -> FldElt
						       
{Compute root of degree one polynomial in variable t}
    
    require Degree(pol) eq 1: "Polynomial must have degree one in t";
    
    num := -Coefficient(pol, 0);
    den := Coefficient(pol, 1);    
    return num/den;
    
end intrinsic;


intrinsic MultipleRoot(fac :: SeqEnum, n :: RngIntElt) -> FldElt
								     
{Compute multiple root of exponent exactly n in the polynomial in t with
specified factorization, assuming it exists and is unique; otherwise, throw
error}
    
    candidates := [f[1]: f in fac | Degree(f[1]) eq 1 and f[2] eq n];

    if #candidates eq 0 then error "No multiple root";
    elif #candidates ge 2 then error "Several multiple roots";
    end if;
    f := candidates[1];
    return DegreeOneRoot(f);
    
end intrinsic;


intrinsic MultipleRoot(pol :: RngUPolElt, n :: RngIntElt) -> FldElt
								    
{Same as above, but factorize on the fly}
    
    fac := Factorization(pol);
    return MultipleRoot(fac, n);
    
end intrinsic;


intrinsic SquareEquation(r::RngElt, v::RngElt) -> Bool, RngElt

{Construct a polynomial equation ensuring that r will be a square,
using variable v}

    pol := v^2 - r;
    fac := Factorization(pol);
    
    if #fac eq 2 then return true, [];
    else return false, pol;
    end if;
    
end intrinsic;


intrinsic SplitCubicEquation(pol :: RngUPolElt, L::SeqEnum[RngElt])
	  -> Bool, SeqEnum[RngElt]
			  
{Equation for polynomial of degree 3 in t to have rational roots; L contains
 variables in which the new equations are written}
    
    require Degree(pol) eq 3: "Must be a cubic equation";

    t := Parent(pol).1;    
    F := Factorization(pol);
    lin := [f[1]: f in F | Degree(f[1]) eq 1];
    
    if #lin eq 3 then
	return true, [];
    elif #lin eq 1 then
	g := [f[1]: f in F | Degree(f[1]) eq 2][1];
	r := Discriminant(g);
	return false, [Denominator(r) * L[1]^2 - Numerator(r)];
    else
	b, r := IsSquare(Discriminant(pol));
	target := Coefficient(pol, 3)
		  * (t-L[1]) * (t-L[2]) * (t-L[3]);
	pol, target;
	if b then //Fix sign of discriminant arbitrarily
	    return false,
		   [Numerator(Coefficient(target, 2) - Coefficient(pol, 2)),
		    Numerator(Coefficient(target, 1) - Coefficient(pol, 1)),
		    Numerator(Coefficient(target, 0) - Coefficient(pol, 0)),
		    r - (L[1]-L[2])*(L[2]-L[3])*(L[3]-L[1])];
	else
	    return false,
		   [Numerator(Coefficient(target, i) - Coefficient(pol, i)):
		    i in [0..2]];
	end if;	
    end if;

end intrinsic;
