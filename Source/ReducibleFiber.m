

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
	F`RatComps := Tate(S, Pl, F`Kodaira);
    end if;
    
    return F;

    /* Ideas: - store root lattice; store dual group; store isomorphism of dual
    group with standard abelian group; store accepted automorphisms on one side
    or the other. */

end intrinsic;


intrinsic HasRationalComponents(S :: EllK3, u :: RngUPolElt, ktype :: MonStgElt)
	  -> Bool, SeqEnum[RngElt]
								       
{Decide if the given elliptic K3 has rational fiber components at the
reducible place u; if not, return one (or maybe two, in case D4)
equations that must be satisfied for the components to be rational}

    require Degree(u) eq 1: "Place must be finite of degree 1";
    
    conf := RootLatticeType(ktype);
    l, n := ParseRootLatticeType(conf);
    
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

    elif ktype eq "IV" then	
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


intrinsic Tate(S :: EllK3, u :: RngUPolElt, ktype :: MonStgElt) ->
	  SeqEnum[SeqEnum[RngUPolElt]]
		 
{Compute fiber components: given a place encoded by a uniformizer u,
compute a sequence of [x0,u^k] or [x0,u^k,y0] describing the
following components: x = x0 + s*u^k, and optionally y = y0 + ... for
some coordinate s when such a parameterization exists. In the case of
place at infinity, return same thing for elliptic surface with t
inverted. In case of fiber of type An, components are cyclically
ordered; in case Dn for n>4, close component comes first.}

    require Degree(u) eq 1: "Place must be finite of degree 1";

    conf := RootLatticeType(ktype);
    l, n := ParseRootLatticeType(conf);

    if l eq "E" then return TateE(S, u, n);
    elif l eq "D" then return TateD(S, u, n);
    elif ktype in ["III", "IV"] then return TateA_add(S, u, ktype);
    else return TateA_mult(S, u, n);
    end if;
    
end intrinsic;


intrinsic TateA_mult(S :: EllK3, u :: RngUPolElt, n :: RngIntElt)
	  -> SeqEnum[SeqEnum[RngUPolElt]]
								       
{Tate in multiplicative case An}

    F := S`BaseField;
    t := S`EllParam;
    Left := [];
    Right := [];
    x0 := F!0;

    //Get first branches: RHS must have a double root
    rhs := RHSQuotientAndEval(S, x0, u, 0, 0);
    x0 := MultipleRoot(rhs, 2);
    
    r := ExactQuotient(rhs, (t-x0)^2);
    r := Evaluate(r, x0);
    b, y0 := IsSquare(r);
    if not b then error "Fibers are not rational, not a square:", r; end if;

    Left := Append(Left, [x0, u, y0*u]);
    Right := Append(Right, [x0, u, -y0*u]);
    k := 1;

    while true do
	//Previous fiber was x = x0 + O(u^k).
	//Adjust x0 such that new fibers on left and right are x = x0+u^(k+1)
	try
	    rhs := RHSQuotientAndEval(S, x0, u, k, 2*k);
	catch e
	    //Reached the end of Tate's algorithm with n even.
	    assert 2*k-2 eq n;
	    break;
	end try;
	
	assert Degree(rhs) eq 2 and Coefficient(rhs, 2) eq r;
	//Does it have a double root? If yes, adjust x0 and add two new branches
	F := Factorization(rhs);
	try
	    x1 := MultipleRoot(F, 2);
	catch e
	    //Reached the end of Tate's algorithm with n odd
	    assert 2*k-1 eq n;
	    //Translate x s.t. equation becomes y^2 = x^2 + (degree 0)
	    num := -Coefficient(rhs, 1);
	    den := 2*Coefficient(rhs, 2);
	    x0 +:= (num/den) * u^k;
	    Left := Append(Left, [x0, u^k]);
	    break;
	end try;
	
	x0 +:= x1*u^k;
	Left := Append(Left, [x0, u^(k+1), y0*u^(k+1)]);
	Right := Append(Right, [x0, u^(k+1), -y0*u^(k+1)]);
	k +:= 1;
    end while;
    
    L := Left cat Reverse(Right);
    return L;
    
end intrinsic;


intrinsic TateA_add(S :: EllK3, u :: RngUPolElt, ktype :: MonStgElt)
	  -> SeqEnum[SeqEnum[RngUPolElt]]
									  
{Tate in additive case An}
    
    F := S`BaseField;
    x0 := F!0;
    
    //Find cube root at u=0
    rhs := RHSQuotientAndEval(S, x0, u, 0, 0);
    x0 := MultipleRoot(rhs, 3);
    //Now look at degree 2
    rhs := RHSQuotientAndEval(S, x0, u, 1, 2);
    if ktype eq "III" then
	assert Degree(rhs) eq 1;
	return [[x0, u]];
    elif ktype eq "IV" then
	assert Degree(rhs) eq 0;
	r := Coefficient(rhs, 0);
	b, y0 := IsSquare(r);
	if not b then error "Fibers are not rational, not a square:", r; end if;
	return [[x0, u, y0*u], [x0, u, -y0*u]];
    else
	error "Wrong Kodaira type", ktype;
    end if;
    
end intrinsic;


intrinsic TateD(S :: EllK3, u :: RngUPolElt, n :: RngIntElt)
	  -> SeqEnum[SeqEnum[RngUPolElt]]
								  
{Tate in case Dn}
    
    F := S`BaseField;
    x0 := F!0;
    L := [];

    //Find cube root at u=0
    rhs := RHSQuotientAndEval(S, x0, u, 0, 0);
    x0 := MultipleRoot(rhs, 3);
    
    //At next step, RHS should have a three simple roots in case D4,
    //one simple and one double root otherwise
    rhs := RHSQuotientAndEval(S, x0, u, 1, 3);
    F := Factorization(rhs);
    
    if n eq 4 then
	lin := [f[1]: f in F
		| Degree(f[1]) eq 1];
	if #lin ne 3 then error "Fibers are not rational, not totally split:", rhs; end if;
	rts := [DegreeOneRoot(f): f in lin];
	for i in [1..3] do
	    L := Append(L, [x0+rts[i]*u, u^2]);
	end for;
	return L;
    end if;

    //Now n is 5 or more
    alpha := MultipleRoot(F, 1);
    beta := MultipleRoot(F, 2);
    //Get close branch
    L := [[x0 + alpha*u, u^2]];
    //Now we need to blow up n-4 times to get the x = x0 + beta*u+... branches
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
	    b, y0 := IsSquare(rhs);
	    if not b then error "Far fibers are not rational, not a square:", rhs; end if;
	    L := Append(L, [x0, u^i, y0*u^i]);
	    L := Append(L, [x0, u^i, -y0*u^i]);
	    break;
	end try;
	
	assert Degree(rhs) eq 2;
	F := Factorization(rhs);
	try
	    x1 := MultipleRoot(F, 2);
	catch e
	    //End computation: this is Dn for n even
	    assert n eq 2*i+2;
	    lin := [f[1]: f in F
		    | Degree(f[1]) eq 1 and f[2] eq 1];
	    if #lin ne 2 then error "Far fibers are not rational, no factorization:", rhs; end if;
	    for f in lin do
		L := Append(L, [x0 + DegreeOneRoot(f)*u^i, u^(i+1)]);
	    end for;
	    break;
	end try;
	
	x0 +:= x1*u^i;
	i +:= 1;
    end while;
    return L;
    
end intrinsic;

intrinsic TateE(S :: EllK3, u :: RngUPolElt, n :: RngIntElt)
	  -> SeqEnum[SeqEnum[RngUPolElt]]
								  
{Tate in case En}

    F := S`BaseField;
    x0 := F!0;
    L := [];

    //Find cube root at u=0
    rhs := RHSQuotientAndEval(S, x0, u, 0, 0);
    x0 := MultipleRoot(rhs, 3);
    //Find cube root at order 1
    rhs := RHSQuotientAndEval(S, x0, u, 1, 3);
    x1 := MultipleRoot(rhs, 3);
    x0 +:= x1*u;
    
    if n eq 6 then
	rhs := RHSQuotientAndEval(S, x0, u, 2, 4);
	assert rhs ne 0 and Degree(rhs) eq 0;
	b, y0 := IsSquare(rhs);
	if not b then error "Far fibers are not rational, not a square:", rhs; end if;
	L := Append(L, [x0, u^2, y0*u^2]);
	L := Append(L, [x0, u^2, -y0*u^2]);
    elif n eq 7 then //Find one far component
	rhs := RHSQuotientAndEval(S, x0, u, 2, 5);
	assert Degree(rhs) eq 1;
	L := Append(L, [x0 + DegreeOneRoot(rhs)*u^2, u^3]);
    else assert n eq 8; //E8: nothing to be done
    end if;
    return L;
    
end intrinsic;
