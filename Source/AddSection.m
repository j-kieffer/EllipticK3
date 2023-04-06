

intrinsic AddSection(~S :: EllK3, x :: RngElt, y :: RngElt)

{Let S know about one of its Mordell-Weil sections}

    require RHS(S,x) eq y^2: "Not a section";
    x := BaseFunctionField(S) ! x;
    y := BaseFunctionField(S) ! y;

    n := Rank(Frame(S));
    M := IntersectionMatrix(S, x, y);
    if IsPositiveDefinite(M) then
        ExtendFrameSpace(~S, M);
        v := FrameSpace(S) ! ([0: i in [1..n]] cat [1]);
        AddSection(~S, v, x, y);
    else
        assert Determinant(M) eq 0;
        v := NSVector(S, [Eltseq(M[n+1])[i]: i in [1..n]]);
        AddSection(~S, v, x, y);
    end if;

end intrinsic;


intrinsic ExtendFrameSpace(~S :: EllK3, int_mat :: AlgMatElt)

{Set frame space of S to quadratic space of rank one more, with inner product
matrix [m,*;*,*] where m is current inner product matrix}

    /* Extend bases of lattices by zero */
    n := Rank(Frame(S));
    r := MordellWeilRank(S);
    int_mat := Matrix(Rationals(), int_mat);
    
    M := Matrix(Rationals(), [Eltseq(b) cat [0]: b in Basis(Frame(S))]);
    S`Frame := Lattice(M, int_mat);
    M := Matrix(Rationals(), [Eltseq(b) cat [0]: b in Basis(NeronSeveriRootSpan(S))]);
    S`RootSpan := Lattice(M, int_mat);
    M := Matrix(Rationals(), [Eltseq(b) cat [0]: b in Basis(RootLattice(S))]);
    S`RootLat := Lattice(M, int_mat);
    M := Matrix(Rationals(), MordellWeilRank(S), n+1,
                [Eltseq(b) cat [0]: b in Basis(MordellWeilLattice(S))]);
    S`MWLat := Lattice(M, int_mat);

    /* Get new quotient groups */
    S`TorsGrp, S`TorsMap := quo < NeronSeveriRootSpan(S) | RootLattice(S) >;
    for i:=1 to #TorsionSections(S) do
        v := NSVector(TorsionSections(S)[i]);
        S`TorsSections[i]`Vec := FrameSpace(S) ! (Eltseq(v) cat [0]);
    end for;
    S`MWGrp, S`MWMap := quo < Frame(S) | RootLattice(S) >;
    for i:=1 to MordellWeilRank(S) do
        v := NSVector(MordellWeilSections(S)[i]);
        S`MWSections[i]`Vec := FrameSpace(S) ! (Eltseq(v) cat [0]);
    end for;

    /* Get new orthogonal projection parallel to RootLat */
    _, mat := OrthogonalizeGram(GramMatrix(RootLattice(S)));
    M := mat * Matrix(Basis(RootLattice(S)));
    n := Rank(RootLattice(S));
    orth_basis := [Frame(S) ! M[i]: i in [1..n]];
    d := Dimension(FrameSpace(S));
    rows := [];
    for i := 1 to d do
        e := [0: k in [1..d]];
        e[i] := 1;
        e := FrameSpace(S) ! e;
        row := e;
        for k := 1 to n do
            v := FrameSpace(S) ! orth_basis[k];
            row := row - InnerProduct(e, v)/InnerProduct(v, v)*v;
        end for;
        Append(~rows, Eltseq(row));
    end for;
    S`MWProj := Matrix(rows);

end intrinsic;


intrinsic NSVector(S :: EllK3, int :: SeqEnum) -> ModTupFldElt

{Return an element in NS(S)\otimes\Q that has the given intersections with
basis elements}

    mat := Matrix(Rationals(), InnerProductMatrix(Frame(S)));
    return Solution(mat, FrameSpace(S)!int);

end intrinsic;


intrinsic NSVector(S :: EllK3, pt :: PtEll) -> ModTupFldElt

{Return an element in NS lattice of S corresponding to the given section}

    n := Rank(Frame(S));
    return Frame(S) ! NSVector(S, Intersections(S, pt)[1..n]);

end intrinsic;


intrinsic AddSection(~S :: EllK3, v :: ModTupFldElt, x :: RngElt, y :: RngElt)

{Let S know about a section (x,y) corresponding to the vector v in its NS ambient space}

    new_frame := ext < Frame(S) | [v] >; // Does the root basis stay the same? Hopefully yes.
    if new_frame eq Frame(S) then // do nothing
        return;
    end if;

    S`Frame := new_frame;
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
            if t eq 0*t then
                continue;
            end if;
            /* Write lift as (element in old NS) + k*v */
            lift := t @@ new_torsmap;
            k := Index(redv, lift @ quomap) - 1;
            assert lift - k*v in NeronSeveriRootSpan(S);
            torssection := New(EllK3MW);
            if lift - k*v in RootLattice(S) then
                torssection`Pt := k*pt;
            else
                pt1 := GenericPoint(S, lift-k*v);
                torssection`Pt := pt1 + k*pt;
            end if;
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
    new_mwgen := FrameSpace(S) ! (v * proj);
    new_mwlat := ext < MordellWeilLattice(S) | [new_mwgen] >;
    new_mwsections := [];
    old_mwbasis := [NSVector(u) * proj: u in MordellWeilSections(S)];
    if Rank(new_mwlat) gt Rank(MordellWeilLattice(S)) then
        /* Can use old basis + v as a new basis */
        new_mwbasis := old_mwbasis cat [v*proj];
        vsection := New(EllK3MW);
        vsection`Pt := pt;
        vsection`Vec := v;
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
                mwpt +:= coefs[i] * GenericPoint(MordellWeilSections(S)[i]);
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


intrinsic AddSection(~S :: EllK3, x :: RngElt)

{Let S know about one of its Mordell-Weil sections}

    rhs := RHS(S, x);
    try
        y := Sqrt(RHS(S, x));
    catch e        
        error "Not the x-coordinate of a section";
    end try;
    AddSection(~S, x, y);
    
end intrinsic;

