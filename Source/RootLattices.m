

intrinsic MonoidOfRootConfigurations() -> GrpAb

{Free abelian group over the symbols in the ADE classification, of rank <= 24}

    names := ["A" cat IntegerToString(i): i in [1..24]]
	     cat ["D" cat IntegerToString(i): i in [4..24]]
	     cat ["E" cat IntegerToString(i): i in [6..8]];
    M := FreeAbelianGroup(48);
    AssignNames(~M, names);
    return M;

end intrinsic;


intrinsic ParseRootLatticeType(conf :: GrpAbElt) -> MonStgElt, RngIntElt

{Return ADE type and rank}

    require Parent(conf) eq MonoidOfRootConfigurations(): "Not a root lattice type";
    S := Eltseq(conf);
    if not (&and[r ge 0: r in S] and &+S eq 1) then
	error "Root configuration must contain a single lattice";
    else
	i := Index(S, 1);
	name := Names(Parent(conf))[i];
	l := name[1];
	n := StringToInteger(name[2..#name]);
	return l, n;
    end if;

end intrinsic;


intrinsic RootLatticeType(name :: MonStgElt) -> GrpAbElt

{Get root lattice type as element of monoid of root configuration}

    // Convert if Kodaira
    name := KodairaToLatticeType(name);
    l := name[1];
    n := StringToInteger(name[2..#name]);
    M := MonoidOfRootConfigurations();
    case l:
    when "A":
	return M.n;
    when "D":
	return M.(n-3+24);
    when "E":
	return M.(n-5+45);
    else
        error "Root lattice name ", name, " not recognized";
    end case;

end intrinsic;


intrinsic RootLattice(name :: MonStgElt) -> Lat

{Standard root lattice of given name An, Dn or En}

    // Convert if Kodaira
    name := KodairaToLatticeType(name);
    l := name[1];
    n := StringToInteger(name[2..#name]);
    
    D := LatticeDatabase();
    L := Lattice(D, name);
    d := Dimension(L);
    G := GramMatrix(L);
    if l eq "D" then // Correct sign of first basis vector
	G[1,3] := -1;
	G[3,1] := -1;
    end if;
    return Lattice(ScalarMatrix(d, 1), G);
    
end intrinsic;


intrinsic RootLattice(type :: GrpAbElt) -> Lat
					      
{Standard root lattice of given type}
    
    S := ElementToSequence(type);
    names := Names(Parent(type));
    L := StandardLattice(0);
    for i := 1 to #names do
	if S[i] lt 0 then error "Root configuration must have positive coefficients";
	end if;
	for j := 1 to S[i] do
	    L := DirectSum(L, RootLattice(names[i]));
	end for;
    end for;
    return L;
    
end intrinsic;


intrinsic RootLatticeCompIndex(type :: GrpAbElt, i :: RngIntElt)
	  -> RngIntElt
		 
{Index of standard basis vector corresponding to simple component number i}

    l, n := ParseRootLatticeType(type);
    if l eq "A" then
	return i;
    elif l eq "D" then
	if i eq 1 then
	    return n;
	elif i eq 2 then
	    return 1;
	elif i eq 3 then
	    return 2;
	end if;
    elif n eq 6 then
	return [1,5][i];
    elif n eq 7 then
	return 6;
    end if;
	
end intrinsic;
