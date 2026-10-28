SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);AttachSpec("mdmagma/v2/mdmagma.spec");tt:=Cputime();F<z>:=GF(2^7);curves:=AssociativeArray();
for j in F do
 for E in Twists(EllipticCurveFromjInvariant(j)) do
  if #E ne 146 then continue;end if;
  mp:=MinimalPolynomial(j,GF(2));assert Degree(mp) eq 7;curves[mp]:=E;
 end for;
end for;assert #Keys(curves) eq 1;
out:="P73/data/quotient19/isogeny3_target_data.m";
PrintFile(out,Sprintf("N:=73;field_poly:=%o;SourceJPolynomials:=[* *];TargetJPolynomials:=[* *];SourceECoefficients:=[* *];",[Integers()!a:a in Eltseq(DefiningPolynomial(F,GF(2)))]):Overwrite:=true);
for mp in Keys(curves) do
 E:=curves[mp];maps:=MDIsogenies(E,3);assert #maps eq 4 and &and[Degree(phi) eq 3:phi in maps];js:=[jInvariant(Codomain(phi)):phi in maps];assert #Seqset(js) eq 4;
 ps:={MinimalPolynomial(j,GF(2)):j in js};assert #ps eq 1;target:=Rep(ps);assert Degree(target) eq 28;
 PrintFile(out,Sprintf("Append(~SourceJPolynomials,%o);Append(~TargetJPolynomials,%o);Append(~SourceECoefficients,%o);",[Integers()!a:a in Eltseq(mp)],[Integers()!a:a in Eltseq(target)],[[Integers()!a:a in Eltseq(c,GF(2))]:c in aInvariants(E)]));
 printf "NONIC73_NATIVE3_TARGET source_j_degree=7 isogenies=4 distinct_targets=4 one_Frob7_cycle=true target_closed_degree=28 norm_multiplicity=1 time=%o\n",Cputime(tt);
end for;
assert (GF(73)!5)^6 eq 3;
PrintFile(out,"print \"NONIC73_NATIVE3_TARGET_DATA_READY\";");print "NONIC73_NATIVE3_TARGETS_COMPLETE";quit;
