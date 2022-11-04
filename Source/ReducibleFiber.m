

intrinsic ReducibleFiber(S :: EllK3, Pl :: RngUPolElt) -> EllK3RedFib

{Reducible fiber of S at the given place}

    if (Pl eq 0) then
	return ReducibleFiber(InvertT(S), S`EllParam);
    end if;

    F := New(EllK3RedFib);
    F`Pl := Pl;
    F`Kodaira := KodairaType(S, Pl);
    F`RootType := RootLatticeType(F`Kodaira);
    
    if Degree(Pl) eq 1 and HasRationalComponents(S, Pl, F`Kodaira) then
	
    end if;

    
    return F;

    /* Ideas: - store root lattice; store dual group; store isomorphism of dual
    group with standard abelian group; store accepted automorphisms on one side
    or the other. */

end intrinsic;


intrinsic HasRationalComponents(S :: EllK3, u :: RngUPolElt, K :: MonStgElt)
	  -> Bool, SeqEnum[RngElt]
								       
{Decide if the given elliptic K3 has rational fiber components at the
reducible place u; if not, return one (or maybe two, in case D4)
equations that must be satisfied for the components to be rational}

    require Degree(u) eq 1: "Place must be finite of degree 1";
    
    conf := RootLatticeType(K);
    l, n := ParseElementaryRootConfiguration(conf);
    
    t := S`EllParam;
    R := Parent(t);
    x0 := R!0;

    B := S`BaseField;
    Q<v1,v2,v3> := PolynomialRing(B, 3);
    L := [v1,v2,v3];
    
    if l eq "E" and n eq 6 then	
	//Find cube root at u=0
	rhs := RHSQuotientAndEval(S, x0, u, 0, 0);
	x0 := MultipleRoot(rhs, 3);
	//Find cube root at order 1
	rhs := RHSQuotientAndEval(S, x0, u, 1, 3);
	x1 := MultipleRoot(rhs, 3);
	x0 +:= x1*u;	
	rhs := RHSQuotientAndEval(S, x0, u, 2, 4);
	assert rhs ne 0 and Degree(rhs) eq 0;
	rhs := Coefficient(rhs, 0);
	return SquareEquation(rhs, L[1]);

    elif l eq "D" and n eq 4 then	
	//Find cube root at u=0
	rhs := RHSQuotientAndEval(S, x0, u, 0, 0);
	x0 := MultipleRoot(rhs, 3);
	//At next step, RHS should have a three simple roots in case D4
	rhs := RHSQuotientAndEval(S, x0, u, 1, 3);
	assert Degree(rhs) eq 3;
	return SplitCubicEquation(rhs, L);

    elif l eq "D" then
	//Find cube root at u=0
	rhs := RHSQuotientAndEval(S, x0, u, 0, 0);
	x0 := MultipleRoot(rhs, 3);
	//At next step, RHS should have one simple and one double root
	rhs := RHSQuotientAndEval(S, x0, u, 1, 3);
	beta := MultipleRoot(rhs, 2);
	x0 +:= beta*u;
	i := 2;
	while true do
	    try
		rhs := RHSQuotientAndEval(S, x0, u, i, 2*i+1);
	    catch e
		//End computation: this is Dn for n odd
		assert n eq 2*i+1;
		rhs := RHSQuotientAndEval(S, x0, u, i, 2*i);
		assert Degree(rhs) eq 0;
		rhs := Coefficient(rhs, 0);
		return SquareEquation(rhs, L[1]);
	    end try;
	    assert Degree(rhs) eq 2;
	    try
		x1 := MultipleRoot(rhs, 2);
	    catch e
		//End computation: this is Dn for n even
		assert n eq 2*i+2;
		D := Discriminant(rhs);
		return SquareEquation(D, L[1]);
	    end try;
	    x0 +:= x1 * u^i;
	    i +:= 1;
	end while;	

    elif K eq "IV" then	
	//Find cube root at u=0
	rhs := RHSQuotientAndEval(S, x0, u, 0, 0);
	x0 := MultipleRoot(rhs, 3);
	//Now look at degree 2
	rhs := RHSQuotientAndEval(S, x0, u, 1, 2);	
	assert Degree(rhs) eq 0;
	r := Coefficient(rhs, 0);
	return SquareEquation(r, L[1]);
	
    elif l eq "A" and n ge 2 then	
	//Get first branches: RHS must have a double root
	rhs := RHSQuotientAndEval(S, x0, u, 0, 0);
	x0 := MultipleRoot(rhs, 2);	
	r := ExactQuotient(rhs, (t-x0)^2);
	r := Evaluate(r, x0);
	return SquareEquation(r, L[1]);
	
    else //Other cases have no automorphisms
	return true, [];
    end if;
end intrinsic;
		       
