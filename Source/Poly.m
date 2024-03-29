

intrinsic Variables(A :: Rng) -> SeqEnum[RngElt]
{Return the list of variables of A}

    return [A.j: j in [1..NumberOfGenerators(A)]];

end intrinsic;


intrinsic DegreeOneRoot(pol :: RngUPolElt) -> FldElt
						       
{Compute root of degree one polynomial in variable t}
    
    require Degree(pol) eq 1: "Polynomial must have degree one in t";
    
    num := -Coefficient(pol, 0);
    den := Coefficient(pol, 1);    
    return num/den;
    
end intrinsic;


intrinsic MultipleRoot(fac :: SeqEnum, n :: RngIntElt) -> FldElt
								     
{Compute multiple root of exponent exactly n in the polynomial in t with
specified factorization, assuming it exists and is unique; otherwise, throw
error}
    
    candidates := [f[1]: f in fac | Degree(f[1]) eq 1 and f[2] eq n];

    if #candidates eq 0 then error "No multiple root";
    elif #candidates ge 2 then error "Several multiple roots";
    end if;
    f := candidates[1];
    return DegreeOneRoot(f);
    
end intrinsic;


intrinsic MultipleRoot(pol :: RngUPolElt, n :: RngIntElt) -> FldElt
								    
{Same as above, but factorize on the fly}
    
    fac := Factorization(pol);
    return MultipleRoot(fac, n);
    
end intrinsic;


intrinsic SquareEquation(r::RngElt, v::RngElt) -> Bool, SeqEnum[RngElt]

{If r is a square, output true and a square root. If not, construct a
polynomial equation ensuring that r will be a square, using variable v}
    
    pol := v^2 - r;
    fac := Factorization(pol);
    
    if #fac eq 2 then return true, DegreeOneRoot(fac[1][1]);
    else return false, [pol];
    end if;
    
end intrinsic;


intrinsic SplitCubicEquation(pol :: RngUPolElt, L::SeqEnum[RngElt])
	  -> Bool, SeqEnum[RngElt]
			  
{Equation for polynomial of degree 3 in t to have rational roots; L contains
 variables in which the new equations are written}
    
    require Degree(pol) eq 3: "Must be a cubic equation";

    t := Parent(pol).1;    
    F := Factorization(pol);
    lin := [f[1]: f in F | Degree(f[1]) eq 1];
    
    if #lin eq 3 then
	return true, [];
    elif #lin eq 1 then
	g := [f[1]: f in F | Degree(f[1]) eq 2][1];
	r := Discriminant(g);
	return false, [Denominator(r) * L[1]^2 - Numerator(r)];
    else
	b, r := IsSquare(Discriminant(pol));
	target := Coefficient(pol, 3)
		  * (t-L[1]) * (t-L[2]) * (t-L[3]);
	pol, target;
	if b then //Fix sign of discriminant arbitrarily
	    return false,
		   [Numerator(Coefficient(target, 2) - Coefficient(pol, 2)),
		    Numerator(Coefficient(target, 1) - Coefficient(pol, 1)),
		    Numerator(Coefficient(target, 0) - Coefficient(pol, 0)),
		    r - (L[1]-L[2])*(L[2]-L[3])*(L[3]-L[1])];
	else
	    return false,
		   [Numerator(Coefficient(target, i) - Coefficient(pol, i)):
		    i in [0..2]];
	end if;	
    end if;

end intrinsic;


intrinsic Reverse(f :: RngUPolElt, n :: RngIntElt) -> RngUPolElt

{Return x^n * f(1/x). The integer n must be at least deg(f).}

    Coeffs := Coefficients(f);
    if #Coeffs gt n+1 then error "Reverse is not a polynomial"; end if;
    Coeffs := Coeffs cat [0: i in [#Coeffs..n]];
    return Parent(f)! Reverse(Coeffs);

end intrinsic;


intrinsic Reverse(f :: FldFunRatUElt, n :: RngIntElt) -> RngUPolElt
{Return x^n * f(1/x)}
    
    t := Parent(f).1;
    return t^n * Evaluate(f, 1/t);

end intrinsic;

intrinsic ReverseCoefficients(f :: RngMPolElt) -> FldFunRatElt

{Given a multivariate polynomial f over k[t], substitute t by 1/t in f}

    C, M := CoefficientsAndMonomials(f);
    deg := Max([0] cat [Degree(c): c in C]);
    res := 0;
    for i:=1 to #C do
        res := Reverse(C[i], deg) * M[i];
    end for;
    t := Parent(f) ! BaseRing(Parent(f)).1;
    return res / t^deg;

end intrinsic;

intrinsic ReverseCoefficients(f :: FldFunRatMElt) -> FldFunRatMElt
                                                                         
{Given a rational fraction f over a base field k[t], substitute t by 1/t in f}

    return ReverseCoefficients(Numerator(f)) / ReverseCoefficients(Denominator(f));

end intrinsic;


intrinsic SeriesExpansion(f :: RngElt, Pl :: RngUPolElt, n :: RngIntElt)
          -> RngUPolElt
                 
{Given a rational function whose denominator doesn't vanish at the given place,
compute its expansion up to O(Pl^n)}

    if Denominator(f) eq 1 then
        return Numerator(f);
    elif Gcd(Denominator(f), Pl) eq 1 then
        _, a, _ := Xgcd(Denominator(f), Pl^n);
        return ((Numerator(f) mod Pl^n) * a) mod Pl^n;
    else
        error "No pole allowed";
    end if;
    
end intrinsic;


intrinsic ReduceCoefficients(f :: RngMPolElt, g :: RngUPolElt) -> RngElt
                                                                      
{Given f, a multivariate polynomial over a univariate polynomial ring, return
the polynomial obtained by reducing all the coefficients of f modulo g}

    P := Parent(f);
    C, M := CoefficientsAndMonomials(f);
    res := P!0;
    for i:=1 to #C do
        res +:= (C[i] mod g) * M[i];
    end for;
    return res;
    
end intrinsic;


intrinsic CoefficientsInT(f :: RngMPolElt) -> SeqEnum

{Given f, a multivariate polynomial over a univariate polynomial ring in t,
return the sequence of coefficients of f with respect to t}

    P := Parent(f);
    C, M := CoefficientsAndMonomials(f);
    m := Max([0] cat [Degree(c): c in C]);
    res := [P| ];
    for i:=0 to m do
        r := P!0;
        for j:=1 to #C do
            r +:= Coefficient(C[j], i) * M[j];
        end for;
        Append(~res, r);
    end for;
    return res;             
    
end intrinsic;

//Magma refuses to do this for some reason.
intrinsic Evaluate(f :: FldFunRatMElt, n :: RngIntElt, x :: FldFunRatMElt) -> FldFunRatMElt
{Replace the n-th variable by x in f}

    return Evaluate(Numerator(f), n, x) / Evaluate(Denominator(f), n, x);
    
end intrinsic;


intrinsic Evaluate(f :: FldFunRatMElt, n :: RngIntElt, x :: RngMPolElt) -> FldFunRatMElt
{Replace the n-th variable by x in f}

    return Evaluate(Numerator(f), n, x) / Evaluate(Denominator(f), n, x);
    
end intrinsic;


intrinsic InvertT(f :: FldFunRatMElt) -> FldFunRatMElt
{Given a rational fraction over k[t], replace t by 1/t}

    Cnum, Mnum := CoefficientsAndMonomials(Numerator(f));
    Cden, Mden := CoefficientsAndMonomials(Denominator(f));
    dnum := Max([0] cat [Degree(c): c in Cnum]);
    dden := Max([0] cat [Degree(c): c in Cden]);
    Cnum := [Reverse(c, dnum): c in Cnum];
    Cden := [Reverse(c, dden): c in Cden];
    t := Parent(f) ! (BaseRing(Parent(f)).1);
    return t^(dden - dnum) * (&+ [Cnum[i] * Mnum[i]: i in [1..#Cnum]])
           / (&+ [Cden[i] * Mden[i]: i in [1..#Cden]]);

end intrinsic;
                                            
