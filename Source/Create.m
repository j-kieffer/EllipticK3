

intrinsic EllipticK3(Coeffs :: SeqEnum[RngUPolElt]:
                     ComputeReducibleFibers := true)
          -> EllK3

{Define an elliptic K3 surface by the data of its coefficients}

    require #Coeffs eq 2 or #Coeffs eq 3: "Must specify either 2 or 3
    coefficients";

    S := New(EllK3);
    F := BaseRing(Coeffs[1]);
    P := Parent(Coeffs[1]);

    /* Convert base ring to a field if not already the case */
    if not IsField(F) then
        FF := FieldOfFractions(F);
        t := P.1;
        P := PolynomialRing(FF);
        try
            names := [Sprint(F.i): i in [1..NumberOfGenerators(F)]];
            AssignNames(~FF, names);
        catch e;
        end try;
        AssignNames(~P, [Sprint(t)]);
    end if;

    /* Convert coefficients to function field in t */
    F := FieldOfFractions(P);
    AssignNames(~F, [Sprint(P.1)]);
    eccoeffs := [F|0,0,0,0,0];
    if #Coeffs eq 2 then
        eccoeffs[4] := Coeffs[1];
        eccoeffs[5] := Coeffs[2];
    else
        eccoeffs[2] := Coeffs[1];
        eccoeffs[4] := Coeffs[2];
        eccoeffs[5] := Coeffs[3];
    end if;

    S`EC := EllipticCurve(eccoeffs);
    if not ComputeReducibleFibers then
        return S;
    end if;

    /* Get reducible fibers */
    list_fib := [];
    Pls := ReduciblePlaces(S);
    for Pl in Pls do
        fib := ReducibleFiber(S, Pl);
        Append(~list_fib, fib);
    end for;
    S`RedFib := list_fib;

    /* Get lattices */
    L := StandardLattice(0);
    S`RootConf := MonoidOfRootConfigurations()!0;
    for fib in S`RedFib do
        L := DirectSum(L, RootLattice(fib));
        S`RootConf +:= RootType(fib);
    end for;
    n := Dimension(L);
    S`Frame := L;
    S`RootSpan := L;
    S`RootLat := L;

    /* Get trivial groups */
    S`TorsGrp, S`TorsMap := quo < S`RootSpan | S`RootLat >;
    S`TorsSections := [];
    S`MWGrp, S`MWMap := quo < S`Frame | S`RootLat >;
    S`MWLat := sub < Frame(S) | [] >;
    S`MWProj := ZeroMatrix(Rationals(), Dimension(L));
    S`MWSections := [];

    return S;

end intrinsic;


intrinsic EllipticK3(L :: Lat) -> EllK3

{Construct the family of elliptic K3 surfaces polarized by the given
Néron--Severi lattice.}

    error "Not implemented";

end intrinsic;


intrinsic InvertT(S :: EllK3:
                  ComputeReducibleFibers := true) -> EllK3

{Make the change of variables t->1/t}

    deg := [4,8,12];
    Coeffs := [Reverse(Coefficients(S)[i], deg[i]): i in [1..3]];
    return EllipticK3(Coeffs: ComputeReducibleFibers := false);

    if ComputeReducibleFibers then
        //Todo: set everything from what we know about S
        error "Not implemented";
    end if;

end intrinsic;
