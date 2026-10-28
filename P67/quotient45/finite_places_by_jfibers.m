// Root-approved replacement after exact reversible pause of2518446/start185376147.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
load "P67/quotient45/finite_common_quick.m";
function PermPower(pp,ii,n)
 for _ in[1..n]do ii:=pp[ii];end for;return ii;
end function;
load "P67/data/quotient45/jbase_relation_modular_data.m";
cc:=[PF!c:c in JBasePolynomialCoefficients];
root:="P67/data/quotient45/";
eligible:=[];jindices:=[];degrees:=[];zs:=[* *];keys:=[* *];
// Prepared fallback only: use certified j-norm fibers for sources as well.
for si in [1..#SourceJPolynomials]do
 pj:=PF!SourceJPolynomials[si];d:=Degree(pj);F<a>:=ext<GF(2)|pj>;PT<V>:=PolynomialRing(F);hp:=&+[(PT!Eltseq(cc[k]))*a^(k-1):k in[1..#cc]];assert hp ne 0;
 rr:=Roots(hp);bps:={MinimalPolynomial(r[1],GF(2)):r in rr};candidates:=Poles(AF!t);
 for pv in bps do candidates cat:=Zeros(Evaluate(pv,AF!t));end for;
 count:=0;
 for pl in candidates do
  if pl in eligible or Degree(pl)ne d then continue;end if;
  // No j-pole has degree7,8,9: the certified pole degrees are1 and11.
  jv:=Evaluate(jAF,pl);if MinimalPolynomial(jv,GF(2))ne pj then continue;end if;
  Append(~eligible,pl);Append(~degrees,d);Append(~jindices,si);z:=PlaceCoordinates(pl,Z);Append(~zs,z);Append(~keys,JKey(z,jv));count+:=1;
 end for;
 assert count eq 11;
 printf "NONIC67H3_SOURCE_JFIBER_COMPLETE index=%o degree=%o places=%o time=%o\n",si,d,count,Cputime(tt);
end for;assert #eligible eq 55 and #Seqset([k:k in keys])eq 55;assert [#[i:i in jindices|i eq k]:k in [1..5]]eq[11,11,11,11,11];perm:=[];
for k in [1..#eligible]do jv:=Evaluate(jAF,eligible[k]);key:=JKey(SigmaCoordinates(zs[k]),jv);ii:=Index(keys,key);assert ii ne 0;Append(~perm,ii);end for;
assert Sort(perm)eq[1..#perm];for i in[1..#perm]do assert PermPower(perm,i,11)eq i;end for;
seen:={};orbs:=[];for i in[1..#eligible]do if i in seen then continue;end if;orb:=[i];k:=perm[i];while k ne i do Append(~orb,k);k:=perm[k];end while;assert #orb eq 11;Append(~orbs,orb);seen join:=Seqset(orb);end for;
gens:=[];for pl in eligible do aa,bb:=TwoGenerators(pl);Append(~gens,[CodeAF(aa),CodeAF(bb)]);end for;
PrintFile(root cat "source_places_data.m",Sprintf("SourcePlaceGenerators:=%o;SourcePlaceKeys:=%o;SourcePlaceDegrees:=%o;SourceJIndices:=%o;SourceSigmaPermutation:=%o;SourceSigmaOrbits:=%o;",gens,keys,degrees,jindices,perm,orbs):Overwrite:=true);
printf "NONIC67H3_SOURCE_PLACES_COMPLETE places=55 sigma_orbits=5 j_multiplicities=[11,11,11,11,11] both_eigenlines_if_applicable=true time=%o\n",Cputime(tt);
for fi in[1..#TargetJPolynomials]do
 pj:=PF!TargetJPolynomials[fi];d:=Degree(pj);F<a>:=ext<GF(2)|pj>;PT<V>:=PolynomialRing(F);hp:=&+[(PT!Eltseq(cc[k]))*a^(k-1):k in [1..#cc]];assert hp ne 0;rr:=Roots(hp);bps:={MinimalPolynomial(r[1],GF(2)):r in rr};candidates:=Poles(AF!t);
 for pv in bps do candidates cat:=Zeros(Evaluate(pv,AF!t));end for;
 targets:=[];ks:=[* *];
 for pl in candidates do
  if pl in targets or Degree(pl)ne d then continue;end if;jv:=Evaluate(jAF,pl);if MinimalPolynomial(jv,GF(2))ne pj then continue;end if;
  Append(~targets,pl);Append(~ks,JKey(PlaceCoordinates(pl,Z),jv));
 end for;
 expected:=11;assert #targets eq expected and #Seqset([k:k in ks])eq expected;pp:=[];
 for k in[1..#targets]do key:=JKey(SigmaCoordinates(PlaceCoordinates(targets[k],Z)),Evaluate(jAF,targets[k]));ix:=Index(ks,key);assert ix ne 0;Append(~pp,ix);end for;
 assert Sort(pp)eq[1..#pp];for k in[1..#pp]do assert PermPower(pp,k,11)eq k;end for;
 tg:=[];for pl in targets do aa,bb:=TwoGenerators(pl);Append(~tg,[CodeAF(aa),CodeAF(bb)]);end for;
 PrintFile(root cat Sprintf("target_places%o_jbase_data.m",fi),Sprintf("TargetPlaceGenerators:=%o;TargetPlaceKeys:=%o;TargetSigmaPermutation:=%o;TargetJPolynomial:=%o;TargetDegree:=%o;",tg,ks,pp,TargetJPolynomials[fi],d):Overwrite:=true);
 printf "NONIC67H3_TARGET_PLACES_COMPLETE fiber=%o degree=%o places=%o finite_base_values=%o infinity_checked=true all_target_places_included=true time=%o\n",fi,d,#targets,#rr,Cputime(tt);
end for;
print "NONIC67H3_FINITE_PLACES_ALL_COMPLETE";

quit;
