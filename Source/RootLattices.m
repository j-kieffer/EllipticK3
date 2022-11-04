

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

    require Parent(conf) eq MonoidOfRootConfigurations(): "Not a root configuration";
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
