
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
    P := PolynomialRing(S);
    a := ExactQuotient(Norm(v), 4);
    num_degree := 2+a;
    Q := PolynomialRing(P, 2*num_degree);
    AssignNames(~Q, ["x", "y"] cat ["a" cat Sprint(i): i in [1..(2*num_degree - 2)]]);
    x := Q.1;
    y := Q.2;
    den := Q!1;
    
    //Non-torsion sections
    d := Dimension(RootLattice(S));
    r := MordellWeilRank(S);
    for i:=d+1 to d+r do
        if v[i] ne 0 then
            error "Not implemented if non-torsion sections are present";
        end if;
    end for;
    assert v in NeronSeveriRootSpan(S);

    //Torsion sections
    if not v in RootLattice(S) then
        tors := [s : s in TorsionSections(S) | NSVector(s) - v in RootLattice(S)][1];
        error "Not implemented if torsion sections are present";
    end if;
    assert v in RootLattice(S);

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
    rest := 2+a - (&+f);
    assert rest ge 0;
    f[1] +:= rest;

    //Compute denominator
    for i:=1 to n do
        F := ReducibleFibers(S)[i];
        if not IsZero(Place(F)) then
            den *:= Q ! Place(F)^(f[i]);
        end if;
    end for;
    
    //Sections take the form (a(t)x + b(t))/den
    section := Q!0;
    t := EllipticParameter(S);
    for i:=0 to num_degree do
        section +:= Q.(i+3) * t^i;
    end for;
    for i := 0 to num_degree - 4 do
        section +:= Q.(i+3 + num_degree + 1) * t^i * x;
    end for;
    section := 1/den * section;

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
        if IsZero(pl) then
            sec := ReverseCoefficients(section);
            sec := Evaluate(sec, x, x/t^4);
            sec := Evaluate(sec, y, y/t^6);
            pl := t;
        else
            sec := section;
        end if;
        for j:=1 to d do
            comp := Component(F, j);
            x0, m := Explode(comp[1..2]);
            val := Max(0, f[i] + f[i]*FiberDivisor(F)[j] - v[k-1+j]);
            ev := SeriesExpansion(sec * pl^(f[i]), pl, val);
            ev := Evaluate(ev, x, x0 + m*x);
            ev := ReduceCoefficients(ev, pl^val);
            //Now ev must be identically zero.
            coefs := CoefficientsInT(ev);
            for c in coefs do
                eqlist cat:= Coefficients(c, x); //what about y?
            end for;
        end for;
        k +:= d;
    end for;
    assert k eq Dimension(RootLattice(S)) + 1;

    //Build linear system
    nrows := 2*num_degree - 2;
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

    //Solve; we should find a 2-dimensional space of sections
    sols := Basis(NullSpace(mat)); //Should be echelonized; simplify?
    assert #sols eq 2;
    //Pick out any non-constant one
    sec1 := Numerator(section);
    sec2 := Numerator(section);    
    for i:=1 to nrows do
        sec1 := Evaluate(sec1, Q.(i+2), sols[1][i]);
        sec2 := Evaluate(sec2, Q.(i+2), sols[2][i]);
    end for;
    if Degree(sec2) eq 0 then
        sec2 := sec1;
    end if;
    sec1 := den;

    //Convert to nicer field
    R<x,y> := RationalFunctionField(PolynomialRing(S), 2);
    s1 := R!0;
    s2 := R!0;
    c1, m1 := CoefficientsAndMonomials(sec1);
    for i:=1 to #c1 do
        k := Degree(m1[i], Q.1);
        l := Degree(m1[i], Q.2);
        s1 +:= c1[i] * x^k * y^l;
    end for;
    c2, m2 := CoefficientsAndMonomials(sec2);
    for i:=1 to #c2 do
        k := Degree(m2[i], Q.1);
        l := Degree(m2[i], Q.2);
        s2 +:= c2[i] * x^k * y^l;
    end for;
    return s2/s1;    
    
end intrinsic;
