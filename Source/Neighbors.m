
intrinsic EvenOrbits(L :: Lat) -> SeqEnum[LatElt]
                                            
{Compute representatives of L/2L modulo Aut(L)}
    
    G := AutomorphismGroup(L);
    G2 := ChangeRing(G, GF(2));
    orbs := OrbitsOfSpaces(G2, 1);
    res := [L ! (Basis(o[2])[1]): o in orbs];
    res := [v: v in res | Norm(v) mod 4 eq 0];
    return res;
    
end intrinsic;

intrinsic FiberDivisor(F :: EllK3RedFib) -> LatElt
                                                
{Return the projection of the fiber divisor onto the root lattice of F}

    //Cf. multiplicities in Tate's algorithm + Magma's standard bases
    l, n := ParseRootLatticeType(RootType(F));
    L := RootLattice(RootType(F));
    if l eq "A" then
        return L ! [-1: i in [1..n]];
    elif l eq "D" then
        return L ! ([-1, -1] cat [-2: i in [1..(n-3)]] cat [-1]);
    elif l eq "E" and n eq 6 then
        return L ! [-1, -2, -3, -2, -1, -2];
    elif l eq "E" and n eq 7 then
        return L ! [-1, -2, -3, -4, -3, -2, -2];
    else //E8
        return L ! [-2, -3, -4, -6, -5, -4, -3, -2];
    end if;

end intrinsic;

intrinsic IsEffective(S :: EllK3, v :: LatElt, a :: RngIntElt) -> Bool, SeqEnum[RngIntElt]
                                                                      
{Return true iff the divisor v + (2+a)F is effective}

    //First, substract non-torsion
    d := Dimension(Frame(S));
    r := MordellWeilRank(S);
    for i := d-r+1 to d do
        if v[i] lt 0 then return false; end if;
        v[i] := 0;
    end for;    
    assert v in NeronSeveriRootSpan(S);
    
    //Then, substract torsion
    if not v in RootLattice(S) then
        ws := [NSVector(sec): sec in TorsionSections(S) | v - NSVector(S) in RootLattice(S)];
        assert #ws eq 1;
        v := v - ws;
    end if;
    assert v in RootLattice(S);

    //Finally, sum up the number of fiber components needed
    n := #ReducibleFibers(S);
    k := 1;
    f := [];
    for i:=1 to n do
        F := ReducibleFibers(S)[i];
        L := RootLattice(KodairaType(F));
        d := Dimension(L);
        w := L ! [v[j]: j in [k..(k+d-1)]];
        w0 := FiberDivisor(F); //negative entries
        Append(~f, Ceiling(Maximum([0] cat [w[j]/w0[j]: j in [1..d]])));
        k +:= d;
    end for;
    assert k eq Dimension(RootLattice(S)) + 1;
    return (&+f le 2+a), f;
    
end intrinsic;

intrinsic EllipticDivisors(L :: Lat, v0 :: LatElt, a :: RngIntElt) -> SeqEnum[LatElt]
                                                                               
{Return the elliptic divisors F' such that F'.F = 2 and F'.O = a-2}
    
    vs := [v[1]: v in CloseVectors(L, -1/2*v0, a, a)];
    //Quotient by action of automorphism group.
    G := AutomorphismGroup(L);
    Gmod2, Gred := ChangeRing(G, GF(2));
    S := Stabilizer(Gmod2, VectorSpace(Gmod2) ! v0) @@ Gred;
    all := [v0 + 2*v: v in vs];  //| IsEffective(S, v0 + 2*v, a)];
    res := [];
    while all ne [] do
        v := all[1];
        orb := Orbit(S, v);
        assert orb subset all;
        all := [w: w in all | not w in orb];
        Append(~res, v);
    end while;
    return res;
    
end intrinsic;

intrinsic IsSubgraph(S :: EllK3, v :: LatElt, f :: SeqEnum[RngIntElt]) -> Bool
                                                     
{Return true iff v corresponds to a subgraph divisor}

    //First, check non-torsion
    L := Frame(S);
    d := Dimension(L);
    r := MordellWeilRank(S);
    for i := d-r+1 to d do
        w := [0: j in [1..d]];
        w[i] := 1;
        if v[i] ne 0 and InnerProduct(v, L!w) ne 0 then
            return false;
        end if;
    end for;
    
    //Then, check torsion sections?    

    //Finally, check fiber components    
    nb := #ReducibleFibers(S);
    k := 1;
    for i:=1 to nb do
        F := ReducibleFibers(S)[i];
        l, n := ParseRootLatticeType(RootType(F));
        for j:=1 to n do
            w := [0: a in [1..d]];
            w[k + (j-1)] := 1;
            t := InnerProduct(v, L!w);
            if (l eq "A" and (j eq 1 or j eq n))
               or (l eq "D" and j eq n-1)
               or (l eq "E" and n eq 6 and j eq 6)
               or (l eq "E" and n eq 7 and j eq 6)
               or (l eq "E" and n eq 8 and j eq 8) then
                //Connects to the extra simple component, but do nothing?
            end if;
            iscomp := (v[k + (j-1)] - f[i] * FiberDivisor(F)[j]) ne 0;
            if iscomp and t ne 0 then
                return false;
            end if;
        end for;
        k +:= n;
    end for;
    assert k eq d+1;
    return true;

end intrinsic;

intrinsic SubgraphDivisors(S :: EllK3, v0 :: LatElt, a :: RngIntElt) -> SeqEnum[LatElt]

{Return the subgraph divisors among the elliptic divisors}

    divs := EllipticDivisors(S, v0, a);
    res := [];
    for v in divs do
        _, f := IsEffective(S, v, a);
        if &+f eq a+2 and IsSubgraph(S, v, f) then
            Append(~res, v);
        end if;
    end for;
    return res;

end intrinsic;

intrinsic NewEllipticParameter(S :: EllK3, v :: LatElt) -> RngElt

{Given an elliptic divisor of the form v + 2O + (2+a)F, compute a new elliptic
parameter for its associated 2-neighbor step.}

    //Setup
    a := ExactQuotient(Norm(v), 4);
    num_degree := 2+a;
    multO := 2;
    
    //Find MW section, if any
    pt := Zero(GenericFiber(S));
    d := Dimension(RootLattice(S));
    r := MordellWeilRank(S);
    secs := MordellWeilSections(S);
    for i:=d+1 to d+r do
        pt +:= v[i] * GenericPoint(secs[i-d]);
    end for;
    if pt ne Zero(GenericFiber(S)) then
        v := v - NSVector(S, pt);
    end if;
    assert v in NeronSeveriRootSpan(S);
    if not v in RootLattice(S) then
        tors := [s : s in TorsionSections(S) | NSVector(s) - v in RootLattice(S)][1];
        pt +:= GenericPoint(tors);
        v := v - NSVector(tors);
    end if;
    assert v in RootLattice(S);
    has_mw := pt ne Zero(GenericFiber(S));
    if has_mw then
        xs, ys, zs := Explode(Coordinates(pt));
        xs := xs/zs;
        ys := ys/zs;
        multO := multO-1; //is that right?
    end if;

    //Distribute (2+a)F to make v + (2+a)F effective
    n := #ReducibleFibers(S);
    k := 1;
    f := [];
    for i:=1 to n do
        F := ReducibleFibers(S)[i];
        L := RootLattice(KodairaType(F));
        d := Dimension(L);
        w := L ! [v[j]: j in [k..(k+d-1)]];
        w0 := FiberDivisor(F); //negative entries
        Append(~f, Ceiling(Maximum([0] cat [w[j]/w0[j]: j in [1..d]])));
        k +:= d;
    end for;
    assert k eq Dimension(RootLattice(S)) + 1;
    rest := num_degree - (&+f);
    assert rest ge 0;
    f[1] +:= rest;

    //Compute weighted monomials
    weights := &cat[[[a,b]: a in [0..Floor(multO/2)]]: b in [0..Floor(multO/3)]];
    weights := [w: w in weights | 2*w[1]+3*w[2] le multO];
    
    //Compute number of variables ai
    nb_var := 2; //x,t
    for w in weights do
        a, b := Explode(w);
        nb_var +:= Max(0, num_degree - 4*a - 6*b + 1);
        if has_mw then
            nb_var +:= Max(0, num_degree - 4*a - 6*b);
        end if;
    end for;

    //Create multivariate fraction field
    Q := RationalFunctionField(PolynomialRing(S), nb_var);
    names := ["x", "y"] cat ["a" cat Sprint(i): i in [0..(nb_var-3)]];
    AssignNames(~Q, names);
    x := Q.1;
    y := Q.2;
    
    //Compute denominator
    denom := PolynomialRing(S)!1;
    for i:=1 to n do
        F := ReducibleFibers(S)[i];
        if not IsZero(Place(F)) then
            denom *:= Place(F)^(f[i]);
        end if;
    end for;
    
    //Compute possible sections in terms of coefficients ai
    section := Q!0;
    t := EllipticParameter(S);
    j := 3;
    for w in weights do
        a, b := Explode(w);
        for i:=0 to num_degree - 4*a - 6*b do
            section +:= Q.j * t^i * x^a * y^b;
            j +:= 1;
        end for;
        if has_mw then
            for i := 0 to num_degree - 4*a - 6*b - 1 do
                section := Q.j * t^i * x^a * y^b * (x-xs)/(y-ys);
                j +:= 1;
            end for;
        end if;
    end for;
    section := section/denom;
    
    //Substract 1 to get a dimension 1 subspace
    d := Degree(denom);
    section := Evaluate(section, 3 + d, 0);

    //Collect linear equations in ai
    eqlist := [];
    k := 1;
    for i:=1 to n do
        F := ReducibleFibers(S)[i];
        pl := Place(F);
        _, d := ParseRootLatticeType(RootType(F));
        if f[i] eq 0 then
            k +:= d;
            continue;
        end if;
        sec := section;
        den := denom;
        if IsZero(pl) then
            sec := InvertT(section);
            sec := Evaluate(sec, 1, Q.1/t^4);
            sec := Evaluate(sec, 2, Q.1/t^6);
            den := Reverse(denom, num_degree);
            pl := t;
        end if;
        for j:=1 to d do
            comp := Component(F, j);
            x0, mx, y0, y1, my := Explode(comp);
            val := Max(0, f[i] + f[i]*FiberDivisor(F)[j] - v[k-1+j]);
            ev := SeriesExpansion(1/den * pl^(f[i]), pl, val) * Numerator(sec);
            xpol := Parent(ev).1;
            ypol := Parent(ev).2;
            ev := Evaluate(ev, 1, x0 + mx*xpol);
            ev := Evaluate(ev, 2, y0 + y1*(x0+mx*xpol) + my*ypol);
            ev := ReduceCoefficients(ev, pl^val);
            //Now ev must be identically zero.
            coefs := CoefficientsInT(ev);
            for c in coefs do
                for d in Coefficients(c, ypol) do
                    eqlist cat:= Coefficients(d, xpol);
                end for;
            end for;
        end for;
        k +:= d;
    end for;
    assert k eq Dimension(RootLattice(S)) + 1;

    //Build linear system
    nrows := nb_var - 2;
    ncols := #eqlist;
    mat := ZeroMatrix(BaseField(S), nrows, ncols);
    for i:=1 to ncols do
        coefs, monomials := CoefficientsAndMonomials(eqlist[i]);
        for j:=1 to #coefs do
            assert Degree(monomials[j]) eq 1;
            k := Index([Q.(l+2) : l in [1..nrows]], monomials[j]);
            mat[k,i] := coefs[j];
        end for;
    end for;
    //Remove zero row
    mat := Matrix(BaseField(S), nrows-1, ncols,
                  [Rows(mat)[i]: i in [1..nrows] | i ne Degree(denom)+1]);

    //Solve; we find a 1-dimensional space of sections
    sols := Basis(NullSpace(mat));
    assert #sols eq 1;
    sol := Eltseq(sols[1]);
    sol := sol[1..Degree(denom)] cat [0] cat sol[Degree(denom)+1..nrows-1];
    sec := section;
    for i:=1 to nrows do
        sec := Evaluate(sec, i+2, sol[i]);
    end for;

    //Convert to nicer field
    R<x,y> := RationalFunctionField(PolynomialRing(S), 2);
    num := R!0;
    den := R!0;
    c, m := CoefficientsAndMonomials(Numerator(sec));
    for i:=1 to #c do
        k := Degree(m[i], xpol);
        l := Degree(m[i], ypol);
        num +:= c[i] * x^k * y^l;
    end for;    
    c, m := CoefficientsAndMonomials(Denominator(sec));
    for i:=1 to #c do
        k := Degree(m[i], xpol);
        l := Degree(m[i], ypol);
        den +:= c[i] * x^k * y^l;
    end for;
    return num/den;
    
end intrinsic;
