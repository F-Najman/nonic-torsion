SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);AttachSpec("mdmagma/v2/mdmagma.spec");tt:=Cputime();
out:="P71/data/quotient26/isogeny3_target_data.m";
PrintFile(out,"N:=71;SourceDegrees:=[];SourceJPolynomials:=[* *];TargetJPolynomials:=[* *];SourceTargetIndices:=[* *];SourceECoefficients:=[* *];SourceFieldPolynomials:=[* *];":Overwrite:=true);
count:=0;targetcount:=0;
for d in [7,8]do
 F<z>:=GF(2^d);curves:=AssociativeArray();
 for j in F do
  if Degree(MinimalPolynomial(j,GF(2))) ne d then continue;end if;
  for E in Twists(EllipticCurveFromjInvariant(j))do
   if #E mod 71 eq 0 then curves[MinimalPolynomial(j,GF(2))]:=E;end if;
  end for;
 end for;
 assert #Keys(curves) eq 1;
 for mp in Keys(curves)do
  E:=curves[mp];maps:=MDIsogenies(E,3);assert #maps eq 4 and &and[Degree(phi) eq 3:phi in maps];js:=[jInvariant(Codomain(phi)):phi in maps];assert #Seqset(js) eq 4;
  ps:=Setseq({MinimalPolynomial(j,GF(2)):j in js});Sort(~ps);assert Sort([Degree(p):p in ps]) eq (d eq 8 select [16,16] else [28]);
  inds:=[];
  for p in ps do
   occ:=#[j:j in js|MinimalPolynomial(j,GF(2)) eq p];assert occ*d eq Degree(p);targetcount+:=1;Append(~inds,targetcount);
   PrintFile(out,Sprintf("Append(~TargetJPolynomials,%o);",[Integers()!a:a in Eltseq(p)]));
  end for;
  PrintFile(out,Sprintf("Append(~SourceDegrees,%o);Append(~SourceJPolynomials,%o);Append(~SourceTargetIndices,%o);Append(~SourceECoefficients,%o);Append(~SourceFieldPolynomials,%o);",d,[Integers()!a:a in Eltseq(mp)],inds,[[Integers()!a:a in Eltseq(c,GF(2))]:c in aInvariants(E)],[Integers()!a:a in Eltseq(DefiningPolynomial(F,GF(2)))]));
  count+:=1;printf "NONIC71_NATIVE3_TARGET source_index=%o degree=%o curve_order=%o maps=4 target_degrees=%o norm_weights_all1=true time=%o\n",count,d,#E,[Degree(p):p in ps],Cputime(tt);
 end for;
end for;assert count eq 2 and targetcount eq 3;
assert (GF(71)!7)^26 eq 3;
PrintFile(out,"print \"NONIC71_NATIVE3_TARGET_DATA_READY\";");print "NONIC71_NATIVE3_TARGETS_COMPLETE";quit;
