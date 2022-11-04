

intrinsic Variables(A :: Rng) -> SeqEnum[RngElt]
{Return the list of variables of A}

    return [A.j: j in [1..NumberOfGenerators(A)]];

end intrinsic;
