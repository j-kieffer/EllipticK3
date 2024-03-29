

intrinsic Intersection(fib :: EllK3RedFib, x :: RngElt, y :: RngElt)
	      -> RngIntElt

{Return index of component where (x,y) intersects the given reducible fiber}

    if Place(fib) eq 0 then
        x := Reverse(x, 4);
        y := Reverse(y, 6);
        Pl := Numerator(Parent(x).1);
    else
        Pl := Place(fib);
    end if;

    /* If infinite, return 0 */
    if IsDivisibleBy(Denominator(x), Pl) then
        return 0;
    end if;
    
    /* Convert to polynomials */
    l, n := ParseRootLatticeType(RootType(fib));
    x := SeriesExpansion(x, Pl, n+3);
    y := SeriesExpansion(y, Pl, n+3); //Add some margin
    
    /* Get correct index in list of fiber components */
    i0 := 0;
    m := 0;
    for i := 1 to #Components(fib) do
	    comp := Component(fib, i);
	    x0, mx, y0, y1, my := Explode(comp);
	    if IsDivisibleBy(x-x0, mx) and IsDivisibleBy(y - y0 - y1*x, my)
           and (mx+my) gt m then
            i0 := i;
            m := (mx+my);
        end if;
    end for;
    return i0;

end intrinsic;

intrinsic IntersectionWithO(S :: EllK3, x :: PtEll) -> RngIntElt
{Compute intersection number of x with the zero section}
    x, _, z := Explode(Eltseq(x));
    return Degree(Denominator(x/z))/2;
end intrinsic;

intrinsic Intersection(S :: EllK3, x1 :: PtEll, x2 :: PtEll) -> RngIntElt
								                                    
{Compute intersection pairing of sections x1 and x2 of S in essential NS
lattice}

    if x1 eq x2 then
	    x1x2 := -2;
    else
        x1x2 := IntersectionWithO(S, x1-x2);
    end if;
    x1_O := IntersectionWithO(S, x1);
    x2_O := IntersectionWithO(S, x2);
    return - 2 - x1_O - x2_O + x1x2;    
    
end intrinsic;


intrinsic Intersections(S :: EllK3, x :: RngElt, y :: RngElt)
	      -> SeqEnum[RngIntElt]

{Return intersection numbers of (x,y) with Néron--Severi basis of S}

    res := [];
    pt := GenericFiber(S) ! [x,y];

    // Get intersections with fiber components
    for fib in ReducibleFibers(S) do
	    i := Intersection(fib, x, y);
	    l, n := ParseRootLatticeType(RootType(fib));
	    int := [0: j in [1..n]];
        if i ge 1 then
            int[i] := -1;
        end if;
	    res cat:= int;
    end for;

    // Get intersections with Mordell-Weil
    for sec in MordellWeilSections(S) do
	    Append(~res, -Intersection(S, pt, GenericPoint(sec)));
    end for;
    Append(~res, 4 + Degree(Denominator(x)));
    
    return res;

end intrinsic;

intrinsic Intersections(S :: EllK3, pt :: PtEll) -> SeqEnum[RngIntElt]

{Return intersection numbers of the given section with Néron--Severi basis of S}
    x, y, z := Explode(Coordinates(pt));
    X := PolynomialRing(S) ! (x/z);
    Y := PolynomialRing(S) ! (y/z);
    return Intersections(S, X, Y);

end intrinsic;

intrinsic IntersectionMatrix(S :: EllK3, x :: RngElt, y :: RngElt)
          -> AlgMatElt
{Return the intersection numbers of (x,y) with a frame basis of S in matrix form}

    x := BaseFunctionField(S)!x;
    y := BaseFunctionField(S)!y;
    
    ints := Intersections(S, x, y);
    n := Dimension(FrameSpace(S));
    old_rows := Rows(InnerProductMatrix(Frame(S)));
    new_rows := [Eltseq(old_rows[i]) cat [ints[i]]: i in [1..n]];
    Append(~new_rows, ints);
    return Matrix(Integers(), n+1, n+1, new_rows);

end intrinsic;
