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
function Reconstruct(r)
 d,i,fp,ecs,pcs:=Explode(r);R<z>:=PolynomialRing(GF(2));F<a>:=ext<GF(2)|R!fp>;
 assert Degree(F) eq d;
 fromvec:=func<v|&+[F!v[k]*a^(k-1):k in [1..#v]]>;
 E:=EllipticCurve([fromvec(v):v in ecs]);P:=E![fromvec(v):v in pcs];
 assert Order(P) eq N;return E,P;
end function;
for d in Ds do
 F<z>:=GF(2^d);keys:={};qualifying:=0;
 for j in F do
  for E in Twists(EllipticCurveFromjInvariant(j)) do
   if #E mod N ne 0 then continue;end if;
   G,g:=AbelianGroup(E);inv:=Invariants(G);m:=inv[#inv];if m mod N ne 0 then continue;end if;
   qualifying+:=1;PN:=(m div N)*g(G.#inv);assert Order(PN) eq N;
   for k in [1..N-1] do
    if GCD(k,N) ne 1 then continue;end if;
    key:=XYKey(E,k*PN);if key[1] eq d then Include(~keys,key);end if;
   end for;
  end for;
 end for;
 union:={};reps:=0;
 for rec in Inputs do
  if rec[1] ne d then continue;end if;
  E,P:=Reconstruct(rec);orb:={XYKey(E,k*P):k in [1..N-1]|GCD(k,N) eq 1};
  assert orb subset keys and IsEmpty(orb meet union);union join:=orb;reps+:=1;
 end for;
 assert union eq keys;
 printf "NONIC_INDEPENDENT_FINITE_COVERAGE N=%o degree=%o all_closed_points=%o source_reps=%o qualifying_twists=%o time=%o\n",N,d,#keys,reps,qualifying,Cputime(tt);
end for;
for rec in Inputs do
 d:=rec[1];i:=rec[2];E,P:=Reconstruct(rec);inputkey:=XYKey(E,P);assert inputkey[1] eq d;
 phs:=MDIsogenies(E,3);assert #phs eq 4;
 native:=AssociativeArray();fast:=AssociativeArray();native_a:=AssociativeArray();fast_a:=AssociativeArray();
 for ph in phs do
  EC:=Codomain(ph);PC:=ph(P);assert Order(PC) eq N;k:=XYKey(EC,PC);
  if not IsDefined(native,k) then native[k]:=0;end if;native[k]+:=1;
  if a ne 1 then ka:=XYKey(EC,a*PC);if not IsDefined(native_a,ka) then native_a[ka]:=0;end if;native_a[ka]+:=1;end if;
 end for;
 ins:=MDIsogeniesFast(E,P,3);assert &+[z[3]:z in ins] eq 4;
 for z in ins do
  EC,PC,m:=Explode(z);assert Order(PC) eq N;k:=XYKey(EC,PC);
  if not IsDefined(fast,k) then fast[k]:=0;end if;fast[k]+:=m;
  if a ne 1 then ka:=XYKey(EC,a*PC);if not IsDefined(fast_a,ka) then fast_a[ka]:=0;end if;fast_a[ka]+:=m;end if;
 end for;
 assert Keys(native) eq Keys(fast);
 for k in Keys(native) do assert native[k] eq fast[k] and (d*native[k]) mod k[1] eq 0;end for;
 if a ne 1 then assert Keys(native_a) eq Keys(fast_a);for k in Keys(native_a) do assert native_a[k] eq fast_a[k] and (d*native_a[k]) mod k[1] eq 0;end for;end if;
 printf "NONIC_INDEPENDENT_HECKE_CONTROL N=%o degree=%o index=%o native_maps=4 fast_orbits=%o weighted_image_degrees=%o observations=%o weights=%o alpha=%o agrees=true time=%o\n",N,d,i,#ins,[k[1]:k in Keys(native)],[native[k]:k in Keys(native)],[d*native[k] div k[1]:k in Keys(native)],a,Cputime(tt);
end for;
printf "NONIC_FINITE_AND_HECKE_CONTROLS_COMPLETE N=%o inputs=%o time=%o\n",N,#Inputs,Cputime(tt);quit;
