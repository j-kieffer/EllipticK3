
intrinsic All2Neighbors(L :: Lat) -> SeqEnum[LatElt]
                                            
{Compute representatives of L/2L modulo Aut(L)}
    
    G := AutomorphismGroup(L);
    G2 := ChangeRing(G, GF(2));
    orbs := OrbitsOfSpaces(G2, 1);
    return [L!0] cat [L ! (Basis(o[2])[1]): o in orbs];
    
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

intrinsic EllipticDivisors(S :: EllK3, v0 :: LatElt, a :: RngIntElt) -> SeqEnum[LatElt]
                                                                               
{Return the elliptic divisors F' such that F'.F = 2 and F'.O = a-2}
    
    vs := [v[1]: v in CloseVectors(Frame(S), -1/2*v0, a, a)];
    return [v0 + 2*v: v in vs | IsEffective(S, v0 + 2*v, a)];
    
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
