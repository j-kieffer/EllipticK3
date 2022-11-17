

intrinsic Intersection(fib :: EllK3RedFib, x :: RngUPolElt, y :: RngUPolElt)
	  -> RngIntElt

{Return index of component where (x,y) intersects the given reducible fiber}

    pl := fib`Pl;

    if pl eq 0 then
	error "Not implemented for fiber at infinity";
    end if;

    if not fib`IsRat then
	error "Not implemented for non-rational fibers";
    end if;
    
    /* Get correct index in list of fiber components */
    i0 := 0;
    for i := 1 to #fib`RatComps do
	comp := fib`RatComps[i];
	x0, mx := Explode(comp[1..2]);
	b, q := IsDivisibleBy(x-x0, mx);
	if b and #comp gt 2 then
	    y0, y1, my := Explode(comp[3..5]);
	    b := IsDivisibleBy(y - y0 - y1*q, my);
	end if;
	if b then
	    i0 := i;
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
    for fib in S`RedFib do
	i := Intersection(fib, x, y);
	l, n := ParseRootLatticeType(fib`RootType);
	int := [0]*n;
	int[i] := 1;
	res cat:= int;
    end for;

    // Get intersections with Mordell-Weil
    for sec in S`MWSections do
	xs, ys, zs := Explode(Eltseq(sec`Pt));
	Append(~res, Intersection(S, xs/zs, ys/zs, x, y));
    end for;

    return res;

end intrinsic;
