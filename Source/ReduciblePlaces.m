

intrinsic ReduciblePlaces(S :: EllK3) -> SeqEnum[RngUPolElt]

 {Return list of reducible places of S; 0 encodes infinity}
     
     t := S`EllParam;
     R := Parent(t);
     D := Discriminant(S);
     
     L := [R|];
     if Degree(D) le 22 then
	 L := Append(L, R!0);
     end if;
     rts := [f[1]: f in Factorization(D) | f[2] ge 2];
     L := L cat rts;
     
     return L;

end intrinsic;
	    

intrinsic KodairaType(S :: EllK3, Pl :: RngUPolElt) -> MonStgElt
							  
{Compute the Kodaira type of fiber of S at the given place}

    if Pl eq 0 then return KodairaType(InvertT(S), S`EllParam); end if;

    D := Discriminant(S);
    c4, c6 := Explode(cInvariants(S));
    vD := Valuation(D, Pl);
    v4 := Valuation(c4, Pl);
    v6 := Valuation(c6, Pl);
    
    if vD eq 0 then return "I0";
    elif vD eq 1 then return "I1";
    elif v4 eq 0 and v6 eq 0 then return "I" cat IntegerToString(vD);
    elif v6 eq 1 then return "II";
    elif v4 eq 1 then return "III";
    elif v6 eq 2 then return "IV";
    elif v4 eq 2 or v6 eq 3 then return "I" cat IntegerToString(vD-6) cat "*";
    elif v6 eq 4 then return "IV*";
    elif v4 eq 3 then return "III*";
    elif v6 eq 5 then return "II*";
    else error "Fiber is too singular, not a K3 surface:", [v4,v6,vD];
    end if;
    
end intrinsic;


intrinsic RootLatticeType(name :: MonStgElt) -> GrpAbElt

{Root lattice contribution of a given Kodaira type}

    /* Convert name to root lattice */
    
    if name eq "II" then error "No lattice contribution";
    elif name eq "III" then name := "A1";
    elif name eq "IV" then name := "A2";
    elif name eq "IV*" then name := "E6";
    elif name eq "III*" then name := "E7";			     
    elif name eq "II*" then name := "E8";
    elif "*" in name then
	assert name[1] eq "I" and name[#name] eq "*";
	n := StringToInteger(name[2..(#name-1)]);	
	name := "D" cat IntegerToString(n+4);
    elif name[1] eq "I" then //I something
	n := StringToInteger(name[2..(#name)]);
	if n lt 2 then error "No lattice contribution";
	else name := "A" cat IntegerToString(n-1);
	end if;
    end if;

    /* Now name is a root lattice type */
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


intrinsic RootLatticeType(S :: EllK3, Pl :: RngUPolElt) -> GrpAbElt

{Root lattice contribution of a given reducible place}
    
    name := KodairaType(S, Pl);
    return RootLatticeType(name);
    
end intrinsic;

