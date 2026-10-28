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
for i in [1..#reps] do
 E:=reps[i][1];P:=reps[i][2];
 printf "NONIC53_MODULI_INPUT index=%o fieldpoly=%o E=%o P=%o xy=%o\n",i,Eltseq(DefiningPolynomial(F,GF(2))),[Eltseq(v,GF(2)):v in aInvariants(E)],[Eltseq(P[k],GF(2)):k in [1..2]],[Eltseq(v,GF(2)):v in TateNormalForm_xy(E,P)];
 phs:=MDIsogenies(E,3);assert #phs eq 4;
 printf "NONIC53_HEC_OUTPUTS index=%o count=%o fields=%o\n",i,#phs,[Degree(BaseRing(Codomain(ph))):ph in phs];
end for;
print "NONIC53_FINITE_KEYS_COMPLETE";quit;
