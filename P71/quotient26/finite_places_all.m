SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
load "P71/quotient26/finite_common.m";
Pj:=MinimalPolynomial(jAF);dj:=Degree(Pj);assert IsDivisibleBy(Degree(AF),dj);dd:=LCM([Denominator(c):c in Eltseq(Pj)]);cc:=[PF!(dd*c):c in Eltseq(Pj)];gg:=GCD(cc);cc:=[ExactQuotient(c,gg):c in cc];dv:=Max([Degree(c):c in cc]);assert dv gt 0 and dv le 360;assert &+[AF!c*jAF^(k-1):k->c in cc] eq 0;
root:="P71/data/quotient26/";
PrintFile(root cat "jbase_relation_data.m",Sprintf("JBaseDegreeJ:=%o;JBaseDegreeParameter:=%o;JBasePolynomialCoefficients:=%o;",dj,dv,[[Integers()!a:a in Eltseq(c)]:c in cc]):Overwrite:=true);
printf "NONIC71_JBASE_RELATION_CERTIFIED degree_j=%o degree_parameter=%o exact_field_identity=true time=%o\n",dj,dv,Cputime(tt);
eligible:=[];jindices:=[];degrees:=[];zs:=[* *];keys:=[* *];
for d in [7,8]do
 pls:=Places(AF,d);
 for pl in pls do
  if Valuation(jAF,pl)lt 0 then continue;end if;jv:=Evaluate(jAF,pl);jj:=[Integers()!c:c in Eltseq(MinimalPolynomial(jv,GF(2)))];ix:=Index([s:s in SourceJPolynomials],jj);if ix eq 0 then continue;end if;
  Append(~eligible,pl);Append(~degrees,d);Append(~jindices,ix);z:=PlaceCoordinates(pl,Z);Append(~zs,z);Append(~keys,JKey(z,jv));
 end for;
end for;assert #eligible eq 10 and #Seqset([k:k in keys])eq 10;assert [#[i:i in jindices|i eq k]:k in [1..2]]eq[5,5];perm:=[];
for k in [1..#eligible]do jv:=Evaluate(jAF,eligible[k]);key:=JKey(SigmaCoordinates(zs[k]),jv);ii:=Index(keys,key);assert ii ne 0;Append(~perm,ii);end for;
assert Sort(perm)eq[1..#perm];for i in[1..#perm]do assert perm[perm[perm[perm[perm[i]]]]]eq i;end for;
seen:={};orbs:=[];for i in[1..#eligible]do if i in seen then continue;end if;orb:=[i];k:=perm[i];while k ne i do Append(~orb,k);k:=perm[k];end while;assert #orb eq 5;Append(~orbs,orb);seen join:=Seqset(orb);end for;
gens:=[];for pl in eligible do aa,bb:=TwoGenerators(pl);Append(~gens,[CodeAF(aa),CodeAF(bb)]);end for;
PrintFile(root cat "source_places_data.m",Sprintf("SourcePlaceGenerators:=%o;SourcePlaceKeys:=%o;SourcePlaceDegrees:=%o;SourceJIndices:=%o;SourceSigmaPermutation:=%o;SourceSigmaOrbits:=%o;",gens,keys,degrees,jindices,perm,orbs):Overwrite:=true);
printf "NONIC71_SOURCE_PLACES_COMPLETE places=10 sigma_orbits=2 j_multiplicities=[5,5] both_eigenlines_if_applicable=true time=%o\n",Cputime(tt);
for fi in[1..#TargetJPolynomials]do
 pj:=PF!TargetJPolynomials[fi];d:=Degree(pj);F<a>:=ext<GF(2)|pj>;PT<V>:=PolynomialRing(F);hp:=&+[(PT!Eltseq(cc[k]))*a^(k-1):k in [1..#cc]];assert hp ne 0;rr:=Roots(hp);bps:={MinimalPolynomial(r[1],GF(2)):r in rr};candidates:=Poles(AF!t);
 for pv in bps do candidates cat:=Zeros(Evaluate(pv,AF!t));end for;
 targets:=[];ks:=[* *];
 for pl in candidates do
  if pl in targets or Degree(pl)ne d or Valuation(jAF,pl)lt 0 then continue;end if;jv:=Evaluate(jAF,pl);if MinimalPolynomial(jv,GF(2))ne pj then continue;end if;
  Append(~targets,pl);Append(~ks,JKey(PlaceCoordinates(pl,Z),jv));
 end for;
 expected:=5;assert #targets eq expected and #Seqset([k:k in ks])eq expected;pp:=[];
 for k in[1..#targets]do key:=JKey(SigmaCoordinates(PlaceCoordinates(targets[k],Z)),Evaluate(jAF,targets[k]));ix:=Index(ks,key);assert ix ne 0;Append(~pp,ix);end for;
 assert Sort(pp)eq[1..#pp];for k in[1..#pp]do assert pp[pp[pp[pp[pp[k]]]]]eq k;end for;
 tg:=[];for pl in targets do aa,bb:=TwoGenerators(pl);Append(~tg,[CodeAF(aa),CodeAF(bb)]);end for;
 PrintFile(root cat Sprintf("target_places%o_jbase_data.m",fi),Sprintf("TargetPlaceGenerators:=%o;TargetPlaceKeys:=%o;TargetSigmaPermutation:=%o;TargetJPolynomial:=%o;TargetDegree:=%o;",tg,ks,pp,TargetJPolynomials[fi],d):Overwrite:=true);
 printf "NONIC71_TARGET_PLACES_COMPLETE fiber=%o degree=%o places=%o finite_base_values=%o infinity_checked=true all_target_places_included=true time=%o\n",fi,d,#targets,#rr,Cputime(tt);
end for;
print "NONIC71_FINITE_PLACES_ALL_COMPLETE";quit;
