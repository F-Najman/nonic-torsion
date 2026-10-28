SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
load "C87/finite_common.m";
load "C87/data/jbase_relation_modular_data.m";
assert #Z eq 33 and Degree(AF)eq 31 and #SourceJPolynomials eq 4;
cc:=[PF!c:c in JBasePolynomialCoefficients];root:="C87/data/";
function FiberPlaces(pj,dd)
 assert Degree(pj)eq dd;F<a>:=ext<GF(2)|pj>;PT<V>:=PolynomialRing(F);hp:=&+[(PT!Eltseq(cc[k]))*a^(k-1):k in[1..#cc]];assert hp ne 0;
 rr:=Roots(hp);bps:={MinimalPolynomial(r[1],GF(2)):r in rr};candidates:=Poles(AF!t);
 for pv in bps do candidates cat:=Zeros(Evaluate(pv,AF!t));end for;
 pls:=[];keys:=[* *];
 for pl in candidates do
  if pl in pls or Degree(pl)ne dd then continue;end if;
  jv:=Evaluate(jAF,pl);if MinimalPolynomial(jv,GF(2))ne pj then continue;end if;
  Append(~pls,pl);Append(~keys,JKey(PlaceCoordinates(pl,Z),jv));
 end for;assert #Seqset([k:k in keys])eq #keys;
 printf "NONIC87_FINITE_FIBER degree=%o places=%o finite_parameter_roots=%o infinity_checked=true time=%o\n",dd,#pls,#rr,Cputime(tt);return pls,keys;
end function;
function CodePlaces(pls)
 result:=[];for pl in pls do aa,bb:=TwoGenerators(pl);Append(~result,[CodeAF(aa),CodeAF(bb)]);end for;return result;
end function;
sourcePlaces:=[* *];sourceKeys:=[* *];sourceGens:=[* *];sourceRep:=[];
for idx in[1..4]do
 pp,keys:=FiberPlaces(PF!SourceJPolynomials[idx],SourceDegrees[idx]);assert #pp eq 4;
 perm:=[];for k in[1..4]do key:=JKey(SigmaCoordinates(PlaceCoordinates(pp[k],Z)),Evaluate(jAF,pp[k]));ii:=Index(keys,key);assert ii ne 0;Append(~perm,ii);end for;
 assert Sort(perm)eq[1..4];o:=[1];k:=perm[1];while k ne 1 do assert k notin o;Append(~o,k);k:=perm[k];end while;assert #o eq 4;
 Append(~sourcePlaces,pp);Append(~sourceKeys,keys);Append(~sourceGens,CodePlaces(pp));Append(~sourceRep,1);
 printf "NONIC87_SOURCE_FIBER_COMPLETE index=%o degree9=true four_diamond_places=true time=%o\n",idx,Cputime(tt);
end for;
PrintFile(root cat "source_places_data.m",Sprintf("SourcePlaceGenerators:=%o;SourcePlaceKeys:=%o;SourceRepresentativeIndices:=%o;",sourceGens,sourceKeys,sourceRep):Overwrite:=true);
targetPlaces:=[* *];targetKeys:=[* *];targetGens:=[* *];
for idx in[1..#TargetJPolynomials]do
 pp,keys:=FiberPlaces(PF!TargetJPolynomials[idx],TargetDegrees[idx]);assert #pp eq 16;
 Append(~targetPlaces,pp);Append(~targetKeys,keys);Append(~targetGens,CodePlaces(pp));
 printf "NONIC87_TARGET_FIBER_COMPLETE index=%o degree54=true all16_levels=true time=%o\n",idx,Cputime(tt);
end for;
PrintFile(root cat "target_places_data.m",Sprintf("TargetPlaceGenerators:=%o;TargetPlaceKeys:=%o;",targetGens,targetKeys):Overwrite:=true);
out:=Open(root cat "direct5_RR_cases.jsonseq","w");survivors:=[];n:=0;
for sid in[1..4]do
 assert #SourceTargetIndices[sid]eq 1;fi:=SourceTargetIndices[sid][1];P:=sourcePlaces[sid][1];assert Degree(P)eq 9;
 for ti in[1..#targetPlaces[fi]]do
  Q:=targetPlaces[fi][ti];assert Degree(Q)eq 54;B:=Q-6*P;assert Degree(B)eq 0;
  printf "NONIC87_DIRECT_A5_RR_START source=%o target=%o time=%o\n",sid,ti,Cputime(tt);dm:=Dimension(B);n+:=1;
  fprintf out,"{\"source_index\":%o,\"target_fiber\":%o,\"target_index\":%o,\"dimension\":%o}\n",sid,fi,ti,dm;Flush(out);
  if dm ne 0 then Append(~survivors,<sid,fi,ti,dm>);end if;
  printf "NONIC87_DIRECT_A5_RR_DONE source=%o target=%o dimension=%o time=%o\n",sid,ti,dm,Cputime(tt);
 end for;
end for;delete out;assert n eq 64;
PrintFile(root cat "direct5_RR_complete.m",Sprintf("Cases:=64;Survivors:=%o;Complete:=true;",survivors):Overwrite:=true);
printf "NONIC87_DIRECT_A5_COMPLETE cases=64 survivors=%o time=%o\n",survivors,Cputime(tt);quit;
