
A<r,s> := PolynomialRing(Rationals(),2);
P<t> := PolynomialRing(A);
a := t*((2*r+1)*t+r);
b := 2*r*s*t^4*(t+1);
c := r*s^2*t^7;
S := EllipticK3([a,b,c]);
RootConfiguration(S);
L := Frame(S);
L;

fp := L![0,0,0,0,0,0,0,-4,-4,-7,-6,-5,-4,-3,-2,-1];
Q, Qmap := quo<L|2*L>;
v0 := (fp@Qmap)@@Qmap;
vp := L!((fp-v0)/2);
assert fp eq v0 + 2*vp;
a := 2;           
assert Norm(fp) eq 4*a;

G := AutomorphismGroup(L);
G2, Gred := ChangeRing(G, GF(2));
orbs := OrbitsOfSpaces(G2,1);
norms := [Norm(L ! (Basis(o[2])[1])): o in orbs];

AttachSpec("spec");
L := DirectSum(RootLattice("E7"), RootLattice("D9"));
v0s := EvenOrbits(L);
#v0s;
for v0 in v0s do print(v0); print(EllipticDivisors(L, v0, 2)); end for;

    
    
    
