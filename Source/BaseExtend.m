
intrinsic BaseExtend(S :: EllK3, K :: Rng:
                     ComputeReducibleFibers:=false) -> EllK3
{Return the base-change of S to the field K}

    require IsField(K): "The new base field must be a field";
    
    P := PolynomialRing(S);
    PP := PolynomialRing(K);
    AssignNames(~PP, Names(P));
    t := PP.1;
    newcoefs := [];
    for c in Coefficients(S) do
        Append(~newcoefs, PP ! [K!x: x in Coefficients(c)]);
    end for;
    S2 := EllipticK3(newcoefs: ComputeReducibleFibers := ComputeReducibleFibers);

    if ComputeReducibleFibers then
        error "Not implemented: add sections here";
    end if;
    return S2;
    
end intrinsic;
