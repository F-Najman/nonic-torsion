SetColumns(0);SetQuitOnError(true);SetMemoryLimit(32*10^9);tt:=Cputime();
load "P73/quotient19/finite_common.m";
load "P73/data/quotient19/jbase_relation_data.m";
assert assigned FiberIndex;i:=StringToInteger(FiberIndex);assert i in [1..1];
p28:=PF!TargetJPolynomials[i];F28<a>:=ext<GF(2)|p28>;PT<V>:=PolynomialRing(F28);
h:=&+[(PT!JBasePolynomialCoefficients[k])*a^(k-1):k in [1..#JBasePolynomialCoefficients]];
assert Degree(h) eq JBaseDegreeParameter;rr:=Roots(h);bps:={MinimalPolynomial(r[1],GF(2)):r in rr};
printf "NONIC73_JV_TARGET_BASE fiber=%o degree_parameter=296 rational_v_roots=%o distinct_base_polynomials=%o degrees=%o time=%o\n",i,#rr,#bps,[Degree(p):p in bps],Cputime(tt);
targets:=[];keys:=[* *];uvs:=[* *];
for pv in bps do
 candidates:=Zeros(Evaluate(pv,AF!t));
 for pl in candidates do
  if Degree(pl) ne 28 or Valuation(jAF,pl) lt 0 then continue;end if;
  jv:=Evaluate(jAF,pl);if MinimalPolynomial(jv,GF(2)) ne p28 then continue;end if;
  if pl in targets then continue;end if;
  Append(~targets,pl);z:=PlaceCoordinates(pl,Z);Append(~keys,JKey(z,jv));
  assert Valuation(uAF,pl) ge 0 and Valuation(vAF,pl) ge 0;
  bb:=Matrix(GF(2),28,28,[c:c in Eltseq(jv^k,GF(2)),k in [0..27]]);assert Rank(bb) eq 28;
  uv:=[Eltseq(Solution(bb,Vector(GF(2),Eltseq(Evaluate(w,pl),GF(2))))):w in [uAF,vAF]];
  Append(~uvs,<[Integers()!c:c in uv[1]],[Integers()!c:c in uv[2]]>);
 end for;
 printf "NONIC73_JV_TARGET_BASE_DONE fiber=%o base_degree=%o candidates=%o targets_so_far=%o time=%o\n",i,Degree(pv),#candidates,#targets,Cputime(tt);
end for;
assert #targets eq 8 and #targets eq #Seqset(targets) and #targets eq #Seqset([k:k in keys]);perm:=[];
for k in [1..#targets]do
 pl:=targets[k];jv:=Evaluate(jAF,pl);key:=JKey(SigmaCoordinates(PlaceCoordinates(pl,Z)),jv);ix:=Index(keys,key);assert ix ne 0;Append(~perm,ix);
end for;assert Sort(perm) eq [1..#perm];
for k in [1..#perm]do assert perm[perm[perm[perm[k]]]] eq k;end for;
gens:=[];for pl in targets do aa,bb:=TwoGenerators(pl);Append(~gens,[CodeAF(aa),CodeAF(bb)]);end for;
out:=Sprintf("P73/data/quotient19/target_places%o_jbase_data.m",i);
PrintFile(out,Sprintf("TargetPlaceGenerators:=%o;TargetPlaceKeys:=%o;TargetSigmaPermutation:=%o;TargetSourceJIndex:=%o;TargetPlaneUV:=%o;TargetJPolynomial:=%o;",gens,keys,perm,i,uvs,TargetJPolynomials[i]):Overwrite:=true);
PrintFile(out,"print \"NONIC73_TARGET_PLACES_JV_READY\";");
printf "NONIC73_TARGET_PLACES_JV_COMPLETE fiber=%o places=%o sigma_permutation=%o all_target_places_included=true time=%o\n",i,#targets,perm,Cputime(tt);quit;
