

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


intrinsic KodairaToLatticeType(ktype :: MonStgElt) -> MonStgElt

{Name of root lattice attached to given Kodaira type}

    if ktype eq "II" then error "No lattice contribution";
    elif ktype eq "III" then name := "A1";
    elif ktype eq "IV" then name := "A2";
    elif ktype eq "IV*" then name := "E6";
    elif ktype eq "III*" then name := "E7";			     
    elif ktype eq "II*" then name := "E8";
    elif "*" in ktype then
	assert ktype[1] eq "I" and ktype[#ktype] eq "*";
	n := StringToInteger(ktype[2..(#ktype-1)]);	
	name := "D" cat IntegerToString(n+4);
    elif ktype[1] eq "I" then //I something
	n := StringToInteger(ktype[2..(#ktype)]);
	if n lt 2 then error "No lattice contribution";
	else name := "A" cat IntegerToString(n-1);
	end if;
    else
	name := ktype;
    end if;

    return name;
    
end intrinsic;


intrinsic RootLatticeType(S :: EllK3, Pl :: RngUPolElt) -> GrpAbElt

{Root lattice contribution of a given reducible place}
    
    return RootLatticeType(KodairaToLatticeType(KodairaType(S, Pl)));
    
end intrinsic;

