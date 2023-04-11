
A<r,s> := PolynomialRing(Rationals(),2);
P<t> := PolynomialRing(A);
a := t*((2*r+1)*t+r);
b := 2*r*s*t^4*(t+1);
c := r*s^2*t^7;
S := EllipticK3([a,b,c]);
RootConfiguration(S);
L := Frame(S);
Basis(L);

p1Fp := L![0,0,0,0,0,0,0,-4,-4,-7,-6,-5,-4,-3,-2,-1];
Q, Qmap := quo<L|2*L>;
v0 := (p1Fp@Qmap)@@Qmap;
vp := L!((p1Fp-v0)/2);
assert p1Fp eq v0 + 2*vp;
a := 2;           
assert Norm(p1Fp) mod 4 eq 0;

G := AutomorphismGroup(L);
G2 := ChangeRing(G, GF(2));
orbs := OrbitsOfSpaces(G2,1);
norms := [Norm(L ! (Basis(o[2])[1])): o in orbs];
vs := All2Neighbors(L);
for v in vs do
    v;
    EllipticDivisors(S, v, 1);
    EllipticDivisors(S, v, 2);
end for;

for v in vs do
    v;
    SubgraphDivisors(S, v, 1);
    SubgraphDivisors(S, v, 2);
end for;
v := L ! [0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1];
