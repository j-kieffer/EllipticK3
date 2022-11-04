

intrinsic Reverse(f :: RngUPolElt, n :: RngIntElt) -> RngUPolElt

{Return x^n * f(1/x). The integer n must be at least deg(f).}

    Coeffs := Coefficients(f);
    if #Coeffs gt n+1 then error "Reverse is not a polynomial"; end if;
    Coeffs := Coeffs cat [0: i in [#Coeffs..n]];
    return Parent(f)! Reverse(Coeffs);

end intrinsic;


intrinsic InvertT(S :: EllK3) -> EllK3

{Make the change of variables t->1/t}
        
    F := S`BaseField;
    t := S`EllParam;
    deg := [4,8,12];
    Coeffs := [Reverse(S`Coeffs[i], deg[i]): i in [1..3]];
    return EllipticK3(Coeffs);
    
end intrinsic;
