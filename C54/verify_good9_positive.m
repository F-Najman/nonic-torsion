SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tp54:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
// Direct Tate reconstruction from the exact native X1(27) coordinate formulas.
data:=eval Read("C54/data/good9_pointdata.txt");def,records:=Explode(data);assert #records gt 0;
R:=PolynomialRing(GF(5));f:=R!def;assert IsIrreducible(f);K<a>:=ext<GF(5)|f>;
units:=[k:k in [1..13]|GCD(k,27) eq 1];allkeys:={};weighted:=0;roots1:=0;roots3:=0;
function PairKeyPositive54(xy)
 return Minimum([Eltseq(xy[1]^(5^i),GF(5)) cat Eltseq(xy[2]^(5^i),GF(5)):i in [0..8]]);
end function;
for i in [1..#records] do
 xycoeff,n,r2,jcoeff,nr3,orbitLength:=Explode(records[i]);xy:=[Evaluate(R!v,a):v in xycoeff];jj:=Evaluate(R!jcoeff,a);
 x,y:=Explode(xy);assert x*y*(x*y-1) ne 0;
 r:=(x^2*y-x*y+y-1)/(x*(x*y-1));s:=(x*y-y+1)/(x*y);b:=r*s*(r-1);c:=s*(r-1);
 EP:=EllipticCurve([1-c,-b,-b,K!0,K!0]);assert Discriminant(EP) ne 0;PP:=EP![0,0];
 assert 27*PP eq EP!0 and 9*PP ne EP!0 and jInvariant(EP) eq jj;
 assert Degree(MinimalPolynomial(jj,GF(5))) in [1,3,9];
 assert LCM([Degree(MinimalPolynomial(z,GF(5))):z in xy]) eq 9 and nr3 eq 2;
 assert TateNormalForm_xy(EP,PP) eq xy;
 b2,b4,b6,b8:=Explode(bInvariants(EP));RR<T>:=PolynomialRing(K);assert #Roots(4*T^3+b2*T^2+2*b4*T+b6) eq r2;
 assert #Roots(DivisionPolynomial(EP,3)) eq nr3;
 keys:=[PairKeyPositive54(TateNormalForm_xy(EP,k*PP)):k in units];assert #Seqset(keys) eq orbitLength and orbitLength in [1,3,9];
 key:=Minimum(keys);assert key notin allkeys;Include(~allkeys,key);weighted+:=orbitLength*r2;
 if r2 eq 1 then roots1+:=1;else assert r2 eq 3;roots3+:=1;end if;
 if i mod 500 eq 0 then print "NONIC54_GOOD9_POSITIVE_PROGRESS",i,weighted,Cputime(tp54);end if;
end for;
assert #allkeys eq #records and weighted eq 218513;
printf "NONIC54_GOOD9_POSITIVE_DONE representatives=%o primitive27=true rational2counts=true orbit_lengths_rechecked=true all_orbits_disjoint=true weighted=%o roots1=%o roots3=%o complete=true time=%o\n",#records,weighted,roots1,roots3,Cputime(tp54);
quit;
