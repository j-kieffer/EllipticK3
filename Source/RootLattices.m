

intrinsic MonoidOfRootConfigurations() -> GrpAb

{Free abelian group over the symbols in the ADE classification, of rank <= 24}

    names := ["A" cat IntegerToString(i): i in [1..24]]
	     cat ["D" cat IntegerToString(i): i in [4..24]]
	     cat ["E" cat IntegerToString(i): i in [6..8]];
    M := FreeAbelianGroup(48);
    AssignNames(~M, names);
    return M;

end intrinsic;


intrinsic RootConfiguration(name:: MonStgElt) -> GrpAbElt

{Return configuration as a group element}
    
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
    else:
    error "Root lattice name ", name, " not recognized";
end case;

end intrinsic;


intrinsic ParseElementaryRootConfiguration(name :: MonStgElt) -> MonStgElt, RngIntElt

{Return ADE type and rank}

    l := name[1];
    n := StringToInteger(name[2..#name]);
    return l, n;

end intrinsic;


intrinsic ParseElementaryRootConfiguration(conf :: GrpAbElt) -> MonStgElt, RngIntElt

{Return ADE type and rank}

    require Parent(conf) eq MonoidOfRootConfigurations(): "Not a root configuration";
    S := Eltseq(conf);
    if not (&and[r ge 0: r in S] and &+S eq 1) then
	error "Root configuration must contain a single lattice";
    else
	i := Index(S, 1);
	return ParseElementaryRootConfiguration(Names(Parent(conf))[i]);
    end if;

end intrinsic;
