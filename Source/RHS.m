

intrinsic RHS(S :: EllK3, x :: RngElt) -> RngElt

{Return right hand side of elliptic surface equation, evaluated at x}

    a2,a4,a6 := Explode(S`Coeffs);
    return x^3+a2*x^2+a4*x+a6;

end intrinsic;
    

intrinsic RHSQuotientAndEval(S :: EllK3,
			     x0 :: RngElt,
			     u :: RngUPolElt,
			     r :: RngIntElt,
			     n :: RngIntElt) -> RngElt
								
{Evaluate right hand side of S at x0 + x*u^r, divide by u^n, evaluate at given
 place, and return result as a polynomial in t}

    require Degree(u) eq 1: "Must be a degree 1 finite place";
			 
    val := DegreeOneRoot(u);
    F := S`BaseField;
    t := S`EllParam;
    R := Parent(t);    
    Z<x> := PolynomialRing(R);
    
    rhs := Coefficients(RHS(S, x0 + x*u^r)); //Sequence of elements in R    
    rhs := [ExactQuotient(c, u^n): c in rhs];
    Ev := [Evaluate(c, val): c in rhs];
    return &+[Ev[i+1] * t^i: i in [0..#Ev-1]];
    
end intrinsic;
