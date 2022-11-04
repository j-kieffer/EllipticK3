
intrinsic Print(S :: EllK3)
{Print S}
    printf "Elliptic K3 surface with parameter %o and coefficients\n", S`EllParam;
    printf "a2 = %o\n", S`Coeffs[1];
    printf "a4 = %o\n", S`Coeffs[2];
    printf "a6 = %o\n", S`Coeffs[3];
    printf "over %o", S`BaseField;
end intrinsic;
    
intrinsic Print(F :: EllK3RedFib)
{Print F}
    print("Reducible fiber of ", F`K3);
    print("at t = ", F`Pl);
    print("Kodaira type: ", F`Kodaira);
end intrinsic;

intrinsic Print(G :: EllK3MW)
{Print G}
    print("Mordell-Weil group of ", G`K3);
    print("isomorphic to ", G`Grp);
end intrinsic;




		
