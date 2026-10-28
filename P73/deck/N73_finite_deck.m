SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");Attach("common/fast_hecke.m");
tt:=Cputime();N:=73;
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

for d in [6,7] do
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
 left:=Keys(lookup);fixed:=0;all:=#left;orbits:=0;
 for k in Keys(lookup) do
  md:=lookup[k];if XYKey(md[1],2*md[2]) eq k then fixed+:=1;end if;
 end for;
 while #left gt 0 do
  k:=Rep(left);md:=lookup[k];orb:={XYKey(md[1],a*md[2]):a in [1..N-1]};
  assert orb subset Keys(lookup);left diff:=orb;orbits+:=1;
 end while;
 if d eq 6 then assert all eq 4 and fixed eq 4 and orbits eq 1;else assert all eq 36 and fixed eq 0 and orbits eq 1;end if;
 printf "NONIC73_FINITE_DECK_CONTROL degree=%o closed_places=%o fixed_by_delta2=%o diamond_orbits=%o complete=true time=%o\n",d,all,fixed,orbits,Cputime(tt);
end for;
print "NONIC73_FINITE_DECK_CONTROL_COMPLETE";quit;
