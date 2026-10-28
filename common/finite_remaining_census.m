SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");Attach("common/fast_hecke.m");
tt:=Cputime();
function XYKey(E,P)
 xx,yy:=Explode(TateNormalForm_xy(E,P));
 d:=LCM(Degree(MinimalPolynomial(xx,GF(2))),Degree(MinimalPolynomial(yy,GF(2))));
 found:=false;
 for bits in [[1,0],[0,1],[1,1]] do
  z:=bits[1]*xx+bits[2]*yy;mp:=MinimalPolynomial(z,GF(2));
  if Degree(mp) eq d then bb:=bits;found:=true;break;end if;
 end for;assert found;
 B:=Matrix(GF(2),[Eltseq(z^i,GF(2)):i in [0..d-1]]);assert Rank(B) eq d;
 cs:=[Eltseq(Solution(B,Vector(GF(2),Eltseq(v,GF(2))))):v in [xx,yy]];
 return [d] cat bb cat [Integers()!c:c in Eltseq(mp)] cat [Integers()!c:c in &cat cs];
end function;


for N in [61,67,71,79,83,89] do
 for d in [7,8,9] do
  F<z>:=GF(2^d);lookup:=AssociativeArray();
  for j in F do
   for E in Twists(EllipticCurveFromjInvariant(j)) do
    if #E mod N ne 0 then continue;end if;
    G,g:=AbelianGroup(E);inv:=Invariants(G);m:=inv[#inv];if m mod N ne 0 then continue;end if;
    PN:=(m div N)*g(G.#inv);assert Order(PN) eq N;
    for k in [1..N-1] do
     P:=k*PN;key:=XYKey(E,P);
     if key[1] eq d then lookup[key]:=<E,P>;end if;
    end for;
   end for;
  end for;
  left:=Keys(lookup);fixed:=0;all:=#left;orbsizes:=[];jds:=[];
  for k in Keys(lookup) do md:=lookup[k];if XYKey(md[1],2*md[2]) eq k then fixed+:=1;end if;end for;
  while #left gt 0 do
   k:=Rep(left);md:=lookup[k];orb:={XYKey(md[1],a*md[2]):a in [1..N-1]};assert orb subset Keys(lookup);left diff:=orb;
   Append(~orbsizes,#orb);Append(~jds,Degree(MinimalPolynomial(jInvariant(md[1]),GF(2))));
  end while;
  printf "NONIC_REMAINING_PRIME_CENSUS N=%o degree=%o closed_places=%o diamond_orbit_sizes=%o j_degrees=%o delta2_fixed=%o complete=true time=%o\n",N,d,all,orbsizes,jds,fixed,Cputime(tt);
 end for;
end for;
print "NONIC_REMAINING_PRIME_CENSUS_COMPLETE";quit;
