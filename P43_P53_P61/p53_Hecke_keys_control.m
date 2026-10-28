SetColumns(0);SetQuitOnError(true);AttachSpec("mdmagma/v2/mdmagma.spec");
F<a>:=GF(512);n:=53;J:={};ces:=[**];jj:={};
for j in F do
 for E in Twists(EllipticCurveFromjInvariant(j)) do
  if #E mod n ne 0 then continue;end if;
  assert #E eq 530;
  assert Degree(MinimalPolynomial(j,GF(2))) eq 9;
  Include(~J,j);Append(~ces,E);
 end for;
end for;
left:=J;reps:=[**];
while #left gt 0 do
 j:=Rep(left);orb:={j^(2^k):k in [0..8]};assert #orb eq 9 and orb subset J;
 E:=Rep([E:E in ces|jInvariant(E) eq j]);G,g:=AbelianGroup(E);P:=10*g(G.1);assert Order(P) eq 53;
 Append(~reps,<E,P>);left diff:=orb;
end while;
printf "NONIC53_FINITE_INPUTS all_j=%o all_twists=%o j_orbits=%o source_closed_places=%o delta_orbits=%o\n",#J,#ces,#reps,#J*26 div 9,#reps;

Attach("common/fast_hecke.m");
function XYKey(E,P)
 xx,yy:=Explode(TateNormalForm_xy(E,P));L:=Parent(xx);
 d:=LCM(Degree(MinimalPolynomial(xx,GF(2))),Degree(MinimalPolynomial(yy,GF(2))));
 found:=false;
 for bits in [[1,0],[0,1],[1,1]] do
  z:=bits[1]*xx+bits[2]*yy;mp:=MinimalPolynomial(z,GF(2));
  if Degree(mp) eq d then bb:=bits;found:=true;break;end if;
 end for;
 assert found;
 B:=Matrix(GF(2),[Eltseq(z^i,GF(2)):i in [0..d-1]]);assert Rank(B) eq d;
 cs:=[Eltseq(Solution(B,Vector(GF(2),Eltseq(v,GF(2))))):v in [xx,yy]];
 return [d] cat bb cat [Integers()!c:c in Eltseq(mp)] cat [Integers()!c:c in &cat cs];
end function;
for i in [1..#reps] do
 E:=reps[i][1];P:=reps[i][2];
 native:=MDIsogenies(E,3);assert #native eq 4;
 nk:=[XYKey(Codomain(ph),ph(P)):ph in native];
 dnk:=[XYKey(Codomain(ph),2*ph(P)):ph in native];
 fast:=MDIsogeniesFast(E,P,3);assert #fast eq 1;
 EE,PP,mul:=Explode(fast[1]);assert mul eq 4 and Degree(BaseRing(EE)) eq 36;
 fk:=XYKey(EE,PP);dfk:=XYKey(EE,2*PP);
 assert Seqset(nk) eq {fk} and Seqset(dnk) eq {dfk};
 assert fk[1] eq 36 and dfk[1] eq 36;
 printf "NONIC53_HEC_KEY_CHECK index=%o native_maps=4 ambient_degree=%o fast_orbits=1 orbit_length=4 image_degree=36 weight=1 diamond2_agrees=true\n",i,[Degree(BaseRing(Codomain(ph))):ph in native];
 printf "NONIC53_ACTUAL_POINT_KEYS index=%o input=%o output=%o diamond2_output=%o\n",i,XYKey(E,P),fk,dfk;
end for;
print "NONIC53_HEC_KEY_COMPARISON_COMPLETE";quit;
