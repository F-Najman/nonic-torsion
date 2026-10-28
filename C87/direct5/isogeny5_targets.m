SetColumns(0);SetQuitOnError(true);SetMemoryLimit(5*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
tt:=Cputime();out:="C87/data/direct5/isogeny5_target_data.m";
PrintFile(out,"N:=87;SourceJPolynomials:=[* *];TargetJPolynomials:=[* *];SourceDegrees:=[];TargetDegrees:=[];SourceTargetIndices:=[* *];SourceECoefficients:=[* *];SourceFieldPolynomials:=[* *];SourceTraces:=[];SourcePointCoefficients:=[* *];SourceGroupInvariants:=[* *];TargetCountsWithinSource:=[* *];TargetNormMultiplicities:=[* *];":Overwrite:=true);
alltargets:=[* *];sourcecount:=0;counts:=[];
for d in [8,9] do
 F<z>:=GF(2^d);curves:=AssociativeArray();pointdata:=AssociativeArray();invs:=AssociativeArray();checked:=0;qualifying:=0;
 for j in F do
  for E in Twists(EllipticCurveFromjInvariant(j)) do
   checked+:=1;n:=#E;
   if n mod 87 ne 0 then continue;end if;
   assert n eq (d eq 8 select 261 else 522);
   G,g:=AbelianGroup(E);iv:=Invariants(G);m:=iv[#iv];assert m mod 87 eq 0;
   PP:=(m div 87)*g(G.#iv);assert Order(PP) eq 87;
   mp:=MinimalPolynomial(j,GF(2));assert Degree(mp) eq d;
   curves[mp]:=E;pointdata[mp]:=PP;invs[mp]:=iv;qualifying+:=1;
  end for;
 end for;
 expected:=d eq 8 select 0 else 4;printf "CENSUS_OBSERVED degree=%o orbits=%o qualifying=%o checked=%o\n",d,#Keys(curves),qualifying,checked;assert #Keys(curves) eq expected and qualifying eq d*#Keys(curves);
 Append(~counts,<#Keys(curves),qualifying,checked>);
 printf "NONIC87_FINITE_CENSUS d=%o j_orbits=%o qualifying_j_twists=%o all_j_twists=%o exact87=true complete=true time=%o\n",d,#Keys(curves),qualifying,checked,Cputime(tt);
 for mp in Keys(curves) do
  E:=curves[mp];PP:=pointdata[mp];maps:=MDIsogenies(E,5);
  assert #maps eq 6 and &and[Degree(phi) eq 5:phi in maps];
  js:=[jInvariant(Codomain(phi)):phi in maps];assert #Seqset(js) eq 6;
  // Verify every pushed level remains exact 87 over the actual kernel field.
  for phi in maps do
   L:=BaseRing(Domain(phi));PL:=Domain(phi)![L!PP[1],L!PP[2],1];assert Order(phi(PL)) eq 87;
  end for;
  ps:=Setseq({MinimalPolynomial(j,GF(2)):j in js});
  assert #ps eq (d eq 8 select 2 else 1);
  inds:=[];cnts:=[];norms:=[];
  for tpoly in ps do
   e:=Degree(tpoly);assert e eq (d eq 8 select 24 else 54);
   cnt:=#[j:j in js|Evaluate(tpoly,j) eq 0];assert cnt eq (d eq 8 select 3 else 6);
   assert d*cnt eq e;
   tl:=[Integers()!a:a in Eltseq(tpoly)];ix:=Index(alltargets,tl);
   if ix eq 0 then Append(~alltargets,tl);ix:=#alltargets;PrintFile(out,Sprintf("Append(~TargetJPolynomials,%o);Append(~TargetDegrees,%o);",tl,e));end if;
   Append(~inds,ix);Append(~cnts,cnt);Append(~norms,1);
  end for;
  sourcecount+:=1;
  PrintFile(out,Sprintf("Append(~SourceJPolynomials,%o);Append(~SourceDegrees,%o);Append(~SourceTargetIndices,%o);Append(~SourceECoefficients,%o);Append(~SourceFieldPolynomials,%o);Append(~SourceTraces,%o);Append(~SourcePointCoefficients,%o);Append(~SourceGroupInvariants,%o);Append(~TargetCountsWithinSource,%o);Append(~TargetNormMultiplicities,%o);",[Integers()!a:a in Eltseq(mp)],d,inds,[[Integers()!a:a in Eltseq(c,GF(2))]:c in aInvariants(E)],[Integers()!a:a in Eltseq(DefiningPolynomial(F,GF(2)))],2^d+1-#E,[[Integers()!a:a in Eltseq(PP[i],GF(2))]:i in [1,2]],invs[mp],cnts,norms));
  printf "NONIC87_NATIVE5_TARGET source=%o jdegree=%o trace=%o isogenies=6 distinct_targets=6 target_j_degrees=%o norm_multiplicities=%o exact87push=true time=%o\n",sourcecount,d,2^d+1-#E,[Degree(p):p in ps],norms,Cputime(tt);
 end for;
end for;
assert sourcecount eq 4;
PrintFile(out,"assert #SourceJPolynomials eq 4;assert #SourceDegrees eq 4;assert #SourceTargetIndices eq 4;print \"NONIC87_NATIVE5_TARGET_DATA_READY\";");
printf "NONIC87_NATIVE5_TARGETS_COMPLETE source_orbits=%o distinct_target_j_orbits=%o census=%o time=%o\n",sourcecount,#alltargets,counts,Cputime(tt);quit;
