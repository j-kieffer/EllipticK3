

intrinsic AddSection(~S :: EllK3, x :: RngUPolElt, y :: RngUPolElt)

{Let S know about one of its Mordell-Weil sections}

    require RHS(S,x) eq y^2: "Not a section";

    M := IntersectionMatrix(S, x, y);
    if IsPositiveDefinite(M) then
        ExtendFrameSpace(~S, M);
        v := FrameSpace(S) ! ([0]*n cat [1]);
        AddSection(~S, v, x, y);
    else
        v := NSVector(S, intersections);
        AddSection(~S, v, x, y);
    end if;

end intrinsic;


intrinsic ExtendFrameSpace(~S :: EllK3, int_mat :: AlgMatElt)

{Set frame space of S to quadratic space of rank one more, with inner product
matrix [m,*;*,*] where m is current inner product matrix}

    /* Extend bases of lattices by zero */
    M := Matrix([Eltseq(b) cat [0]: b in Basis(Frame(S))]);
    S`Frame := Lattice(M, int_mat);
    M := Matrix([Eltseq(b) cat [0]: b in Basis(NeronSeveriRootSpan(S))]);
    S`RootSpan := Lattice(M, int_mat);
    M := Matrix([Eltseq(b) cat [0]: b in Basis(RootLattice(S))]);
    S`RootLat := Lattice(M, int_mat);
    M := Matrix([Eltseq(b) cat [0]: b in Basis(MordellWeilLattice(S))]);
    S`MWLat := Lattice(M, int_mat);

    /* Get new quotient groups */
    S`TorsGrp, S`TorsMap := quo < NeronSeveriRootSpan(S) | RootLattice(S) >;
    for i:=1 to #TorsionSections(S) do
        v := NSVector(TorsionSections(S)[i]);
        S`TorsSections[i]`Vec := FrameSpace(S) ! (Eltseq(v) cat [0]);
    end for;
    S`MWGrp, S`MWMap := quo < Frame(S) | RootLattice(S) >;
    for i:=1 to MordellWeilRank(S) do
        v := NSVector(TorsionSections(S)[i]);
        S`MWSections[i]`Vec := FrameSpace(S) ! (Eltseq(v) cat [0]);
    end for;

    /* Get new orthogonal projection parallel to RootLat */
    _, mat := OrthogonalizeGram(GramMatrix(RootLattice(S)));
    orth_basis := [Frame(S) ! v : v in mat * Matrix(Basis(RootLattice(S)))];
    d := Dimension(FrameSpace(S));
    rows := [];
    for i := 1 to d do
        e := [0]*d;
        e[i] := 1;
        e := FrameSpace(S) ! e;
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

    return Frame(S) ! NSVector(S, Intersections(S, pt));

end intrinsic;


intrinsic AddSection(~S :: EllK3, v :: ModTupFldElt, x :: RngUPolElt, y :: RngUPolElt)

{Let S know about a section (x,y) corresponding to the vector v in its NS ambient space}

    new_frame := ext < Frame(S) | [v] >; // Does the root basis stay the same? Hopefully yes.
    if new_frame eq Frame(S) then // do nothing
        return;
    end if;

    pt := EllipticCurve(S) ! [x,y];

    /* Update torsion sections */
    new_rootspan := DualBasisLattice(RootLattice(S)) meet new_frame;
    if new_rootspan ne NeronSeveriRootSpan(S) then
        new_torsgrp, new_torsmap := quo < new_rootspan | RootLattice(S) >;
        
        /* List representatives of quotient new/old as reductions of k*v */
        quogroup, quomap := quo < new_rootspan | NeronSeveriRootSpan(S) >;
        vv := new_rootspan!v;
        ordv := Order(vv@quomap); 
        redv := [k*(vv@quomap): k in [0..(ordv-1)]];        
        new_torssections := [];
        
        for t in new_torsgrp do
            /* Write lift as (element in old NS) + k*v */
            k := Index(redv, (t @@ new_torsmap) @ quomap);
            pt1 := Section(S, lift-k*v); // As point on elliptic curve over function field
            torssection := New(EllK3MW);
            torssection`Pt := pt1 + k*pt;
            torssection`Vec := NSVector(S, torssection`Pt);
            Append(~new_torssections, torssection);
        end for;
        /* Reduce size of representatives? */
        
        S`RootSpan := new_rootspan;
        S`TorsGrp := new_torsgrp;
        S`TorsMap := new_torsmap;
        S`TorsSections := new_torssections;
    end if;

    /* Update MW group */
    new_mwgrp, new_mwmap := quo < new_frame | RootLattice(S) >;
    S`MWGrp := new_mwgrp;
    S`MWMap := new_mwmap;

    /* Update MW lattice and generators */
    proj := MordellWeilProjection(S);
    new_mwgen := v * proj;
    new_mwlat := ext < MordellWeilLattice(S) | [new_mwgen] >;
    new_mwsections := [];
    old_mwbasis := [NSVector(u) * proj: u in MordellWeilSections(S)];
    if Rank(new_mwlat) gt Rank(MordellWeilLattice(S)) then
        /* Can use old basis + v as a new basis */
        new_mwbasis := old_mwbasis cat [v*proj];
        vsection := New(EllK3MW);
        vsection`Pt := pt;
        vsection`vec := v;
        new_mwsections := MordellWeilSections(S) cat [vsection];
        S`MWLat := new_mwlat;
        S`MWSections := new_mwsections;
    elif new_mwlat ne MordellWeilLattice(S) then
        /* Similar to torsion */
        new_mwbasis := Basis(new_mwlat);
        quogroup, quomap := quo < new_mwlat | MordellWeilLattice(S) >;
        vv := new_mwlat ! (v * proj);
        ordv := Order(vv@quomap);
        redv := [k*(vv@quomap): k in [0..(ordv-1)]];        
        for b in new_mwbasis do
            k := Index(redv, b@quomap);
            coefs := Coordinates(MordellWeilLattice(S) ! (b - k*new_mwgen));
            mwpt := Zero(EllipticCurve(S));
            for i := 1 to #old_mwbasis do
                mwpt +:= coefs[i] * Section(MordellWeilSections(S)[i]);
            end for;
            mwpt +:= k*pt;

            mwsection := New(EllK3MW);
            mwsection`Pt := mwpt;
            mwsection`Vec := NSVector(S, mwpt);
            Append(~new_mwsections, mwsection);
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
        error "Not the x-coordinate of a section";
    else
        AddSection(~S, x, y);
    end if;

end intrinsic;

