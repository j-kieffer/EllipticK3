

intrinsic Intersection(S :: EllK3, fib :: EllK3RedFib,
                       x :: RngUPolElt, y :: RngUPolElt)
	      -> RngIntElt

{Return index of component where (x,y) intersects the given reducible fiber}

    require FieldOfRationality(fib) eq BaseField(S): "Not implemented for
    non-rational fibers";

    if Place(fib) eq 0 then
        x := Reverse(x, 4);
        y := Reverse(y, 6);
    end if;
    
    /* Get correct index in list of fiber components */
    i0 := 0;
    m := 0;
    for i := 1 to #Components(fib) do
	    comp := Component(fib, i);
	    x0, mx := Explode(comp[1..2]);
	    b, q := IsDivisibleBy(x-x0, mx);
	    if b and #comp gt 2 then
	        y0, y1, my := Explode(comp[3..5]);
	        b := IsDivisibleBy(y - y0 - y1*x, my);
	    end if;
	    if b and mx gt m then
	        i0 := i;
            m := mx;
	    end if;
    end for;
    return i0;

end intrinsic;


intrinsic Intersection(S :: EllK3, x1 :: PtEll, x2 :: PtEll) -> RngIntElt
								                                    
{Compute intersection pairing of sections x1 and x2 of S in essential NS
lattice}

    if x1 eq x2 then
	    x1x2 := -2;
    else
	    x,_,_ := Explode(Coordinates(x1-x2));
	    int := Degree(Denominator(x));
    end if;
    x1_O := Degree(Denominator(x1));
    x2_O := Degree(Denominator(x2));
    return 2 + x1_O + x2_O - x1x2;    
    
end intrinsic;


intrinsic Intersections(S :: EllK3, x :: RngUPolElt, y :: RngUPolElt)
	      -> SeqEnum[RngIntElt]

{Return intersection numbers of (x,y) with Néron--Severi basis of S}

    res := [];

    // Get intersections with fiber components
    for fib in ReducibleFibers(S) do
	    i := Intersection(S, fib, x, y);
	    l, n := ParseRootLatticeType(RootType(fib));
	    int := [0: i in [1..n]];
        if i gt 0 then
	        int[i] := -1;
        end if;
	    res cat:= int;
    end for;

    // Get intersections with Mordell-Weil
    for sec in MordellWeilSections(S) do
	    xs, ys, zs := Explode(Eltseq(sec`Pt));
	    Append(~res, -Intersection(S, xs/zs, ys/zs, x, y));
    end for;
    Append(~res, Degree(Numerator(x)));
    
    return res;

end intrinsic;


intrinsic IntersectionMatrix(S :: EllK3, x :: RngUPolElt, y :: RngUPolElt)
          -> AlgMatElt
{Return the intersection numbers of (x,y) with a frame basis of S in matrix form}
    
    ints := Intersections(S, x, y);
    n := Dimension(FrameSpace(S));
    old_rows := Rows(InnerProductMatrix(Frame(S)));
    new_rows := [Eltseq(old_rows[i]) cat [ints[i]]: i in [1..n]];
    Append(~new_rows, ints);
    return Matrix(new_rows);

end intrinsic;
