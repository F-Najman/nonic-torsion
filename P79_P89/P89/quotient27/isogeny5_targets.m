SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);AttachSpec("mdmagma/v2/mdmagma.spec");tt:=Cputime();F<z>:=GF(2^9);curves:=AssociativeArray();
for j in F do
 for E in Twists(EllipticCurveFromjInvariant(j)) do
  if #E ne 534 then continue;end if;
  mp:=MinimalPolynomial(j,GF(2));assert Degree(mp) eq 9;curves[mp]:=E;
 end for;
end for;assert #Keys(curves) eq 3;
out:="P79_P89/data/P89/quotient27/isogeny5_target_data.m";
PrintFile(out,Sprintf("N:=89;field_poly:=%o;SourceJPolynomials:=[* *];TargetJPolynomials:=[* *];SourceECoefficients:=[* *];",[Integers()!a:a in Eltseq(DefiningPolynomial(F,GF(2)))]):Overwrite:=true);
for mp in Keys(curves) do
 E:=curves[mp];maps:=MDIsogenies(E,5);assert #maps eq 6 and &and[Degree(phi) eq 5:phi in maps];js:=[jInvariant(Codomain(phi)):phi in maps];assert #Seqset(js) eq 6;
 ps:={MinimalPolynomial(j,GF(2)):j in js};assert #ps eq 1;target:=Rep(ps);assert Degree(target) eq 54;
 PrintFile(out,Sprintf("Append(~SourceJPolynomials,%o);Append(~TargetJPolynomials,%o);Append(~SourceECoefficients,%o);",[Integers()!a:a in Eltseq(mp)],[Integers()!a:a in Eltseq(target)],[[Integers()!a:a in Eltseq(c,GF(2))]:c in aInvariants(E)]));
 printf "NONIC89_NATIVE5_TARGET source_j_degree=9 isogenies=6 distinct_targets=6 one_Frob9_cycle=true target_closed_degree=54 norm_multiplicity=1 time=%o\n",Cputime(tt);
end for;
assert (GF(89)!3)^70 eq 5;
PrintFile(out,"print \"NONIC89_NATIVE5_TARGET_DATA_READY\";");print "NONIC89_NATIVE5_TARGETS_COMPLETE";quit;
