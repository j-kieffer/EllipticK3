

intrinsic ReducibleFiber(S :: EllK3, Pl :: RngUPolElt) -> EllK3RedFib

{Reducible fiber of S at the given place}

    require Pl eq 0 or Degree(Pl) ge 1: "Not a valid place";
    Pl := PolynomialRing(S)!Pl;

    if (Pl eq 0) then
        t := EllipticParameter(S);
	    F := ReducibleFiber(InvertT(S), t);
	    F`Pl := 0*t;
        return F;
    else
        F := New(EllK3RedFib);
        F`Pl := Pl;
        F`Kodaira := KodairaType(S, Pl);
        b, eqs := HasRationalComponents(S, Pl, KodairaType(F));
        if b then
            F`Fld := BaseField(S);
	        F`Comps := Tate(S, Pl, KodairaType(F));
        elif #eqs eq 1 then
            if BaseField(S) eq Rationals() or IsNumberField(BaseField(S)) then
                K := NumberField(eqs[1]);
            else
                K := FieldOfFractions(quo<PolynomialRing(S) | eqs[1]>);
            end if;            
            AssignNames(~K, ["a"]);
            S2 := BaseExtend(S, K: ComputeReducibleFibers:=false);
            F`Fld := K;
            F`Comps := Tate(S2, PolynomialRing(S2)!Pl, KodairaType(F));
        end if;
        //Todo: implement D4 in case of non-rational components.
        return F;
    end if;
    
end intrinsic;


intrinsic RHSQuotientAndEval(S :: EllK3,
			                 x0 :: RngElt,
			                 u :: RngUPolElt,
			                 r :: RngIntElt,
			                 n :: RngIntElt: add := 0)
          -> RngElt
								                    
{Evaluate right hand side of S at x0 + x*u^r, divide by u^n, evaluate at given
 place, and return result as a polynomial in t}

    require Degree(u) eq 1: "Must be a degree 1 finite place";
	
    val := DegreeOneRoot(u);
    t := EllipticParameter(S);
    R := Parent(t);    
    Z<x> := PolynomialRing(R);
    
    rhs := Coefficients(RHS(S, x0 + x*u^r) + add); //Sequence of elements in R
    rhs := [ExactQuotient(c, u^n): c in rhs];
    Ev := [Evaluate(c, val): c in rhs];
    return &+[Ev[i+1] * t^i: i in [0..#Ev-1]];
    
end intrinsic;


intrinsic HasRationalComponents(S :: EllK3, u :: RngUPolElt, ktype :: MonStgElt)
	      -> Bool, SeqEnum[RngElt]
						  
{Decide if the given elliptic K3 has rational fiber components at the
reducible place u; if not, return one (or maybe two, in case D4)
equations that must be satisfied for the components to be rational}

    require Degree(u) eq 1: "Place must be finite of degree 1";

    u := PolynomialRing(S)!u;
    conf := RootLatticeType(ktype);
    l, n := ParseRootLatticeType(conf);
    
    t := EllipticParameter(S);
    R := Parent(t);
    x0 := R!0;

    B := BaseField(S);
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
	    return SquareEquation(rhs, t);

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
		        return SquareEquation(rhs, t);
	        end try;
	        assert Degree(rhs) eq 2;
	        try
		        x1 := MultipleRoot(rhs, 2);
	        catch e
		        //End computation: this is Dn for n even
		        assert n eq 2*i+2;
		        D := Discriminant(rhs);
		        return SquareEquation(D, t);
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
	    return SquareEquation(r, t);
	    
    elif l eq "A" and n ge 2 then	
	    //Get first branches: RHS must have a double root
	    rhs := RHSQuotientAndEval(S, x0, u, 0, 0);
	    x0 := MultipleRoot(rhs, 2);	
	    r := ExactQuotient(rhs, (t-x0)^2);
	    r := Evaluate(r, x0);
	    return SquareEquation(r, t);
	    
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
    L := [];

    if l eq "E" then L := TateE(S, u, n);
    elif l eq "D" then L := TateD(S, u, n);
    elif ktype in ["III", "IV"] then L := TateA_add(S, u, ktype);
    else L := TateA_mult(S, u, n);
    end if;
    assert &and [#c eq 5: c in L];
    assert #L eq n;
    return L;
    
end intrinsic;


intrinsic TateA_mult(S :: EllK3, u :: RngUPolElt, n :: RngIntElt)
	      -> SeqEnum[SeqEnum[RngUPolElt]]
					
{Tate in multiplicative case An}

    F := BaseField(S);
    t := EllipticParameter(S);
    Left := [];
    Right := [];
    x0 := F!0;

    //Get first branches: RHS must have a double root
    rhs := RHSQuotientAndEval(S, x0, u, 0, 0);
    x0 := MultipleRoot(rhs, 2);
    
    r := ExactQuotient(rhs, (t-x0)^2);
    r := Evaluate(r, x0);
    b, y0 := IsSquare(r);
    if not b and n ge 2 then
        error "Fibers are not rational, not a square:", r;
    end if;
    k := 1;
	rhs := RHSQuotientAndEval(S, x0, u, k, 2*k);

    while true do
	    //Previous fiber was x = x0 + O(u^k).
	    //Adjust x0 such that new fibers on left and right are x = x0+u^(k+1)
	    
	    assert Degree(rhs) eq 2 and Coefficient(rhs, 2) eq r;
	    //Does it have a double root? If yes, adjust x0 and add two new branches
	    F := Factorization(rhs);
	    try
	        x1 := MultipleRoot(F, 2);
	    catch e
	        //Reached the end of Tate's algorithm with n odd
	        assert 2*k-1 eq n;
	        //Translate x s.t. equation becomes y^2 = x^2 + (degree 0)
	        //num := -Coefficient(rhs, 1);
	        //den := 2*Coefficient(rhs, 2);
	        //x0 +:= (num/den) * u^k;
	        Left := Append(Left, [x0, u^k, 0, 0, u^k]);
	        break;
	    end try;
        
        try
            rhs := RHSQuotientAndEval(S, x0+x1*u^k, u, k+1, 2*(k+1));
        catch e
            //Reached the end with n even
            assert 2*k eq n;	    
	        Left := Append(Left, [x0, u^k, -y0*(x1*u^k+x0), y0, u^(k+1)]);
	        Right := Append(Right, [x0, u^k, y0*(x1*u^k+x0), -y0, u^(k+1)]);
            break;
        end try;
        
        Left := Append(Left, [x0, u^k, -y0*(x1*u^k+x0), y0, u^(k+1)]);
        Right := Append(Right, [x0, u^k, y0*(x1*u^k+x0), -y0, u^(k+1)]);
        x0 +:= x1*u^k;
	    k +:= 1;
    end while;
    
    L := Left cat Reverse(Right);
    return L;
    
end intrinsic;


intrinsic TateA_add(S :: EllK3, u :: RngUPolElt, ktype :: MonStgElt)
	      -> SeqEnum[SeqEnum[RngUPolElt]]
					
{Tate in additive case An}
    
    F := BaseField(S);
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
	    return [[x0, u, y0*u, 0, u^2], [x0, u, -y0*u, 0, u^2]];
    else
	    error "Wrong Kodaira type", ktype;
    end if;
    
end intrinsic;


intrinsic TateD(S :: EllK3, u :: RngUPolElt, n :: RngIntElt)
	      -> SeqEnum[SeqEnum[RngUPolElt]]
					
{Tate in case Dn}
    
    F := BaseField(S);
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
        L := [[x0 + rts[1]*u, u^2, 0, 0, u^2],
              [x0 + rts[2]*u, u^2, 0, 0, u^2],
              [x0, u, 0, 0, u^2],
              [x0 + rts[3]*u, u^2, 0, 0, u^2]];
	    return L;
    end if;

    //Now n is 5 or more
    alpha := MultipleRoot(F, 1);
    beta := MultipleRoot(F, 2);
    L := [[PolynomialRing(S)!0] : i in [1..n]];
    //Get close branch
    L[n] := [x0 + alpha*u, u^2, 0, 0, u^2];
    L[n-1] := [x0, u, 0, 0, u^2];
    //Now we need to blow up n-4 times to get the x = x0 + beta*u+... branches
    x0 +:= beta*u;
    i := 2;
    
    while true do
        L[n-2*i+2] := [x0, u^i, 0, 0, u^i];
	    try
	        rhs := RHSQuotientAndEval(S, x0, u, i, 2*i+1);
	    catch e
	        //End computation: this is Dn for n odd
	        assert n eq 2*i+1;
	        rhs := RHSQuotientAndEval(S, x0, u, i, 2*i);
	        assert Degree(rhs) eq 0;
	        b, y0 := IsSquare(rhs);
	        if not b then error "Far fibers are not rational, not a square:", rhs; end if;
	        L[1] := [x0, u^i, y0*u^i, 0, u^(i+1)];
	        L[2] := [x0, u^i, -y0*u^i, 0, u^(i+1)];
	        break;
	    end try;
	    assert Degree(rhs) eq 2;

        L[n-2*i+1] := [x0, u^i, 0, 0, u^(i+1)];
	    F := Factorization(rhs);
	    try
	        x1 := MultipleRoot(F, 2);
	    catch e
	        //End computation: this is Dn for n even
	        assert n eq 2*i+2;
	        lin := [f[1]: f in F
		            | Degree(f[1]) eq 1 and f[2] eq 1];
	        if #lin ne 2 then error "Far fibers are not rational, no factorization:", rhs; end if;
            L[1] := [x0 + DegreeOneRoot(lin[1]) * u^i, u^(i+1), 0, 0, u^(i+1)];
            L[2] := [x0 + DegreeOneRoot(lin[2]) * u^i, u^(i+1), 0, 0, u^(i+1)];
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

    F := BaseField(S);
    x0 := F!0;
    L := [];

    //Find cube root at u=0
    rhs := RHSQuotientAndEval(S, x0, u, 0, 0);
    x0 := MultipleRoot(rhs, 3);
    node1 := [x0, u, 0, 0, u^2];
    //Find cube root at order 1
    rhs := RHSQuotientAndEval(S, x0, u, 1, 3);
    x1 := MultipleRoot(rhs, 3);
    x0 +:= x1*u;
    node2 := [x0, u^2, 0, 0, u^2];
    
    if n eq 6 then
	    rhs := RHSQuotientAndEval(S, x0, u, 2, 4);
	    assert rhs ne 0 and Degree(rhs) eq 0;
	    b, y0 := IsSquare(rhs);
	    if not b then error "Far fibers are not rational, not a square:", rhs; end if;
        y0 := BaseField(S) ! y0;
        //Get linear equation for the far fibers
        lin := RHSQuotientAndEval(S, x0, u, 2, 5: add := -y0^2*u^4);
        assert Degree(lin) le 1;
        c := Coefficient(lin, 1);
        d := Coefficient(lin, 0);
        L := [[x0, u^2, y0*u^2 - c/(2*y0)*x0*u + d/(2*y0)*u^3, c/(2*y0)*u, u^4],
              [x0, u^2, y0*u^2, 0, u^3],
              node2,
              [x0, u^2, -y0*u^2, 0, u^3],
              [x0, u^2, -y0*u^2 + c/(2*y0)*x0*u - d/(2*y0)*u^3, -c/(2*y0)*u, u^4],
              node1];
        
    elif n eq 7 then //Fix this; far fiber is a parabola, what happens in the middle?
	    rhs := RHSQuotientAndEval(S, x0, u, 2, 5);
	    assert Degree(rhs) eq 1;
        a := DegreeOneRoot(rhs);
        L := [[x0 + a*u^2, u^3, 0, 0, u^3],
              [x0, u^2, 0, 0, u^3],
              [x0, u^2, 0, 0, u^3], //These two should really involve a: does it matter?
              [x0, u^2, 0, 0, u^3],
              node2,
              node1,
              [x0, u^2, 0, 0, u^3]];
        
    else assert n eq 8;
         L := [[x0, u^2, 0, 0, u^3],
               [x0, u^2, 0, 0, u^3],
               [x0, u^2, 0, 0, u^3],
               [x0, u^2, 0, 0, u^3],
               [x0, u^2, 0, 0, u^3],
               [x0, u^2, 0, 0, u^3],
               node2,
               node1];         
    end if;
    return L;
    
end intrinsic;
