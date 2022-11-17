

intrinsic AddSection(~S :: EllK3, x :: RngUPolElt, y :: RngUPolElt)

{Let S know about one of its Mordell-Weil sections}

    require Parent(x) eq Parent(S`EllParam) and Parent(y) eq Parent(S`EllParam): "Not an element of the correct polynomial ring";
    require RHS(S,x) eq y^2: "Not a section";

    /* Get new intersection matrix */
    ints = Intersections(S, x, y);
    n := Dimension(S`Ambient);
    old_rows := Rows(InnerProductMatrix(S`NSLat));
    new_rows := [Eltseq(r[i]) cat [ints[i]]: i in [1..n]];
    Append(~new_rows, intersections);
    int_mat := Matrix(new_rows);

    if IsPositiveDefinite(new_int_mat) then
	ExtendAmbient(~S, int_mat);
	v := S`Ambient ! ([0]*n cat [1]);
	AddNSVector(~S, v, x, y);
    else
	v := NSVector(S, intersections);
	assert InnerProduct(v, v) eq 2;
	AddNSVector(~S, v, x, y);
    end if;

end intrinsic;


intrinsic ExtendAmbient(~S :: EllK3, int_mat :: AlgMatElt)

{Set ambient space of S to quadratic space with inner product matrix [m,*;*,*]
where m is current inner product matrix}

    /* Extend bases of lattices by zero */
    M := Matrix([Eltseq(b) cat [0]: b in Basis(S`NSLat)]);
    S`NSLat := Lattice(M, int_mat);
    M := Matrix([Eltseq(b) cat [0]: b in Basis(S`RootNSLat)]);
    S`RootNSLat := Lattice(M, int_mat);
    M := Matrix([Eltseq(b) cat [0]: b in Basis(S`RootLat)]);
    S`RootLat := Lattice(M, int_mat);
    M := Matrix([Eltseq(b) cat [0]: b in Basis(S`RootDual)]);
    S`RootDual := Lattice(M, int_mat);
    M := Matrix([Eltseq(b) cat [0]: b in Basis(S`MWLat)]);
    S`MWLat := Lattice(M, int_mat);
    
    S`Ambient := AmbientSpace(S`NSLat);    

    /* Get new quotient groups */
    S`TorsGrp, S`TorsMap := quo < S`RootNSLat | S`RootLat >;
    for i:=1 to #S`TorsSections do
	v := S`TorsSections[i]`Vec;
	S`TorsSections[i]`Vec := S`Ambient ! (Eltseq(v) cat [0]);
    end for;
    S`MWGrp, S`MWMap := quo < S`NSLat | S`RootLat >;
    for i:=1 to #S`MWSections do
	v := S`MWSections[i]`Vec;
	S`MWSections[i]`Vec := S`Ambient ! (Eltseq(v) cat [0]);
    end for;
    
    /* Get new orthogonal projection parallel to RootLat */
    _, mat := OrthogonalizeGram(GramMatrix(S`RootLat));
    orth_basis := [S`Ambient!v : v in mat * Matrix(Basis(S`RootLat))];
    d := Dimension(S`Ambient);
    rows := [];
    for i := 1 to d do
	e := [0]*d;
	e[i] := 1;
	e := S`Ambient!e;
	row := e;
	for v in orth_basis do
	    row := row - InnerProduct(e, v)/InnerProduct(v, v)*v;
	    Append(~rows, Eltseq(row));
	end for;
    end for;
    S`MWProj := Matrix(rows);    

end intrinsic;


intrinsic NSVector(S :: EllK3, int :: SeqEnum[FldRatElt]) -> ModTupFldElt

{Return an element in NS(S)\otimes\Q that has the given intersections with
basis elements}

    mat := Matrix(Rationals(), InnerProductMatrix(L));
    return Solution(mat, S`Ambient!int);

end intrinsic;


intrinsic NSVector(S :: EllK3, pt :: PtEll) -> ModTupFldElt

{Return an element in NS lattice of S corresponding to the given section}

    return S`NSLat ! NSVector(S, Intersections(S, pt));

end intrinsic;


intrinsic AddNSVector(~S :: EllK3, v :: ModTupFldElt, x :: RngUPolElt, y :: RngUPolElt)

{Let S know about a section (x,y) corresponding to the vector v in its NS ambient space}

    new_nslat := ext < S`NSLat | [v] >; // Does the root basis stay the same? Hopefully yes.
    if new_nslat eq S`NSLat then // do nothing
	return;
    end if;
    
    pt := S`EC ! [x,y];

    /* Update torsion sections */
    new_rootnslat := S`RootDual meet new_nslat;    
    if new_rootnslat ne S`RootNSLat then
	new_torsgrp, new_torsmap := quo < new_rootnslat | S`RootLat >;
	new_torssections := [];

	// Do it more efficiently: we know order of v modulo old torsgroup...
	for t in new_torsgrp do
	    lift := t @@ new_torsmap;
	    // Write lift as (element in old NS) + k*v
	    k := 0;
	    while not lift - k*v in S`NSLat do
		k +:= 1;
	    end while;
	    pt1 := Section(S, lift-k*v); // As point on elliptic curve over function field
	    torssection := New(EllK3MW);
	    torssection`Pt := pt1 + k*pt;
	    torssection`Vec := NSVector(S, torssection`Pt);
	    Append(~new_torssections, torssection);
	end for;
	
	S`RootNSLat := new_rootnslat;
	S`TorsGrp := new_torsgrp;
	S`TorsMap := new_torsmap;
	S`TorsSections := new_torssections;
    end if;

    /* Update MW group */
    new_mwgrp, new_mwmap := quo < new_nslat | S`RootLat >;
    S`MWGrp := new_mwgrp;
    S`MWMap := new_mwmap;

    /* Update MW lattice and generators */
    new_mwgen := v * S`MWProj;
    new_mwlat := ext < S`MWLat | [new_mwgen] >;
    if new_mwlat ne S`MWLat then
	new_mwsections := [];
	old_mwbasis := [u`Vec * S`MWProj: u in S`MWSections];
	new_mwbasis := Basis(new_mwlat);

	for b in new_mwbasis do
	    // Write b as (element in old MW) + k*new_mwgen
	    k := 0;
	    while not b - k*new_mwgen in S`MWLat do
		k +:= 1;
	    end while;
	    coefs := Coordinates(S`MWLat ! (b - k*new_mwgen));
	    mwpt := Zero(S`EC);
	    for i := 1 to #old_mwbasis do
		mwpt +:= coefs[i] * S`MWSections[i]`Pt;
	    end for;
	    mwpt +:= k*pt;
	    
	    mwsection := New(EllK3MW);
	    mwsection`Pt := mwpt;
	    mwsection`Vec := NSVector(S, mwsection`Pt);
	end for;

	S`MWLat := new_mwlat;
	S`MWSections := new_mwsections;
    end if;
    
end intrinsic;


intrinsic AddSection(~S :: EllK3, x :: RngUPolElt)

{Let S know about one of its Mordell-Weil sections}

    rhs := RHS(S, x);
    R<v> := PolynomialRing(Parent(x));
    b, y := SquareEquation(rhs, v);
    if not b then
	error "Right hand side is not a square";
    else
	AddSection(~S, x, y);
    end if;

end intrinsic;
    
