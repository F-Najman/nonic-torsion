// Standalone finite-field/model validator and re-emitter. No RR or normalization.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "P67/data/quotient45/canonical_data.m";
load "P67/data/quotient45/jmap_q_data.m";
load "P67/data/quotient45/plane_model_certified_data.m";
load "P67/data/quotient15/isogeny3_target_data.m";
assert Submatrix(qcoeffs1400,1,1,45,320)eq ChangeRing(qcoeffs,GF(2));
PP:=ProjectiveSpace(GF(2),44);PR:=CoordinateRing(PP);mons:=[PR.p[1]*PR.p[2]:p in pairs];eqs:=[&+[row[i]*mons[i]:i in[1..1035]]:row in Rows(QuadRelations)];CC:=Curve(PP,eqs);assert Dimension(CC)eq 1 and Degree(CC)eq 88 and ArithmeticGenus(CC)eq 45;
printf "NONIC67_FRESH_MODEL_BINDING canonical_degree88=true arithmetic_genus45=true quadrics903=true time=%o\n",Cputime(tt);
PF:=PolynomialRing(GF(2));H:={1,66,38,29,37,30};assert #H eq 6;finitecounts:=[* *];
function Normalize(z)
 k:=Index([a ne 0:a in z],true);assert k ne 0;return[a/z[k]:a in z];
end function;
function CodeKey(z,j)
 return <[Integers()!a:a in Eltseq(MinimalPolynomial(j,GF(2)))],[[Integers()!a:a in Eltseq(v,GF(2))]:v in Normalize(z)]>;
end function;
function ValidateKey(key)
 pj:=PF!key[1];assert IsIrreducible(pj);d:=Degree(pj);assert d in[7,8,9,24,28,36]and (2^d mod 67)notin H;F<j>:=ext<GF(2)|pj>;z:=[F!PF!a:a in key[2]];assert #z eq 45 and z eq Normalize(z);
 quadvals:=Vector(F,[z[p[1]]*z[p[2]]:p in pairs]);assert IsZero(quadvals*Transpose(ChangeRing(QuadRelations,F)));
 pows:=[*[F|1]:i in[1..45]*];for i in[1..45]do for k in[1..10]do Append(~pows[i],pows[i][#pows[i]]*z[i]);end for;end for;
 ms:=[&*[pows[i][es[i]+1]:i in[1..45]]:es in Std10Exponents];nu:=&+[F!NumeratorVec[i]*ms[i]:i in[1..#ms]];de:=&+[F!DenominatorVec[i]*ms[i]:i in[1..#ms]];assert de ne 0 and nu eq j*de;
 lp:=[&+[F!ProjectionRows[i,k]*z[k]:k in[1..45]]:i in[1..3]];assert lp ne [F|0,0,0];assert &+[lp[2]^es[1]*lp[3]^es[2]*lp[1]^(plane_degree-es[1]-es[2]):es in PlaneExponents]eq 0;
 sz:=Normalize([&+[F!Diamond2Matrix[i,k]*z[k]:k in[1..45]]:i in[1..45]]);return CodeKey(z,j),CodeKey(sz,j);
end function;
// Bind source and target polynomials to fresh untouched elliptic3-isogenies.
for si in[1..5]do
 F<a>:=ext<GF(2)|PF!SourceFieldPolynomials[si]>;E:=EllipticCurve([F!PF!v:v in SourceECoefficients[si]]);assert Degree(F)eq SourceDegrees[si]and #E mod 67 eq 0 and MinimalPolynomial(jInvariant(E),GF(2))eq PF!SourceJPolynomials[si];ph:=MDIsogenies(E,3);assert #ph eq 4 and &and[Degree(p)eq 3:p in ph];jj:=[jInvariant(Codomain(p)):p in ph];assert #Seqset(jj)eq 4;polys:=Sort(Setseq({MinimalPolynomial(j,GF(2)):j in jj}));assert polys eq Sort([PF!TargetJPolynomials[k]:k in SourceTargetIndices[si]]);
 for pj in polys do assert #[j:j in jj|MinimalPolynomial(j,GF(2))eq pj]*SourceDegrees[si]eq Degree(pj);end for;
end for;
print "NONIC67_FRESH_ELLIPTIC_SOURCE_AND_ALL3_ISOGENY_PATTERNS_COMPLETE";
load "P67/data/quotient45/source_places_data.m";
assert #SourcePlaceKeys eq 55;newkeys:=[* *];sigmakeys:=[* *];
for key in SourcePlaceKeys do z,s:=ValidateKey(key);Append(~newkeys,z);Append(~sigmakeys,s);end for;assert newkeys eq SourcePlaceKeys and #Seqset([k:k in newkeys])eq 55;
perm:=[Index(newkeys,k):k in sigmakeys];assert perm eq SourceSigmaPermutation and Sort(perm)eq[1..55];orbs:=[];seen:={};for k in[1..55]do if k in seen then continue;end if;oo:=[k];m:=perm[k];while m ne k do Append(~oo,m);m:=perm[m];end while;assert #oo eq 11;Append(~orbs,oo);seen join:=Seqset(oo);end for;assert orbs eq SourceSigmaOrbits and #orbs eq 5;
for si in[1..5]do ix:=[i:i in[1..55]|SourceJIndices[i]eq si];assert #ix eq 11 and &and[newkeys[i][1]eq SourceJPolynomials[si]:i in ix];Append(~finitecounts,<SourceDegrees[si],2^SourceDegrees[si]mod 67,11>);end for;
out:="P67/data/quotient45/finite_validation_20261026/validated_keys.m";PrintFile(out,Sprintf("ValidatedSourceKeys:=%o;ValidatedSourcePermutation:=%o;ValidatedSourceOrbits:=%o;ValidatedTargetKeys:=[**];ValidatedTargetPermutations:=[**];",newkeys,perm,orbs):Overwrite:=true);
printf "NONIC67_FRESH_SOURCE_KEYS_COMPLETE count55=true distinct=true source_orbits5_by11=true time=%o\n",Cputime(tt);
load "P67/data/quotient45/target_places1_jbase_data.m";
assert #TargetPlaceKeys eq 11;kk:=[* *];ss:=[* *];for key in TargetPlaceKeys do z,sz:=ValidateKey(key);Append(~kk,z);Append(~ss,sz);end for;assert kk eq TargetPlaceKeys and #Seqset([k:k in kk])eq 11;tp:=[Index(kk,k):k in ss];assert tp eq TargetSigmaPermutation and Sort(tp)eq[1..11];assert &and[key[1]eq TargetJPolynomials[1]:key in kk];Append(~finitecounts,<TargetDegree,2^TargetDegree mod 67,11>);PrintFile(out,Sprintf("Append(~ValidatedTargetKeys,%o);Append(~ValidatedTargetPermutations,%o);",kk,tp));printf "NONIC67_FRESH_TARGET_KEYS_COMPLETE fiber=1 degree=%o count11=true distinct=true Frobenius_line_count11=true time=%o\n",TargetDegree,Cputime(tt);
load "P67/data/quotient45/target_places2_jbase_data.m";
assert #TargetPlaceKeys eq 11;kk:=[* *];ss:=[* *];for key in TargetPlaceKeys do z,sz:=ValidateKey(key);Append(~kk,z);Append(~ss,sz);end for;assert kk eq TargetPlaceKeys and #Seqset([k:k in kk])eq 11;tp:=[Index(kk,k):k in ss];assert tp eq TargetSigmaPermutation and Sort(tp)eq[1..11];assert &and[key[1]eq TargetJPolynomials[2]:key in kk];Append(~finitecounts,<TargetDegree,2^TargetDegree mod 67,11>);PrintFile(out,Sprintf("Append(~ValidatedTargetKeys,%o);Append(~ValidatedTargetPermutations,%o);",kk,tp));printf "NONIC67_FRESH_TARGET_KEYS_COMPLETE fiber=2 degree=%o count11=true distinct=true Frobenius_line_count11=true time=%o\n",TargetDegree,Cputime(tt);
load "P67/data/quotient45/target_places3_jbase_data.m";
assert #TargetPlaceKeys eq 11;kk:=[* *];ss:=[* *];for key in TargetPlaceKeys do z,sz:=ValidateKey(key);Append(~kk,z);Append(~ss,sz);end for;assert kk eq TargetPlaceKeys and #Seqset([k:k in kk])eq 11;tp:=[Index(kk,k):k in ss];assert tp eq TargetSigmaPermutation and Sort(tp)eq[1..11];assert &and[key[1]eq TargetJPolynomials[3]:key in kk];Append(~finitecounts,<TargetDegree,2^TargetDegree mod 67,11>);PrintFile(out,Sprintf("Append(~ValidatedTargetKeys,%o);Append(~ValidatedTargetPermutations,%o);",kk,tp));printf "NONIC67_FRESH_TARGET_KEYS_COMPLETE fiber=3 degree=%o count11=true distinct=true Frobenius_line_count11=true time=%o\n",TargetDegree,Cputime(tt);
load "P67/data/quotient45/target_places4_jbase_data.m";
assert #TargetPlaceKeys eq 11;kk:=[* *];ss:=[* *];for key in TargetPlaceKeys do z,sz:=ValidateKey(key);Append(~kk,z);Append(~ss,sz);end for;assert kk eq TargetPlaceKeys and #Seqset([k:k in kk])eq 11;tp:=[Index(kk,k):k in ss];assert tp eq TargetSigmaPermutation and Sort(tp)eq[1..11];assert &and[key[1]eq TargetJPolynomials[4]:key in kk];Append(~finitecounts,<TargetDegree,2^TargetDegree mod 67,11>);PrintFile(out,Sprintf("Append(~ValidatedTargetKeys,%o);Append(~ValidatedTargetPermutations,%o);",kk,tp));printf "NONIC67_FRESH_TARGET_KEYS_COMPLETE fiber=4 degree=%o count11=true distinct=true Frobenius_line_count11=true time=%o\n",TargetDegree,Cputime(tt);
load "P67/data/quotient45/target_places5_jbase_data.m";
assert #TargetPlaceKeys eq 11;kk:=[* *];ss:=[* *];for key in TargetPlaceKeys do z,sz:=ValidateKey(key);Append(~kk,z);Append(~ss,sz);end for;assert kk eq TargetPlaceKeys and #Seqset([k:k in kk])eq 11;tp:=[Index(kk,k):k in ss];assert tp eq TargetSigmaPermutation and Sort(tp)eq[1..11];assert &and[key[1]eq TargetJPolynomials[5]:key in kk];Append(~finitecounts,<TargetDegree,2^TargetDegree mod 67,11>);PrintFile(out,Sprintf("Append(~ValidatedTargetKeys,%o);Append(~ValidatedTargetPermutations,%o);",kk,tp));printf "NONIC67_FRESH_TARGET_KEYS_COMPLETE fiber=5 degree=%o count11=true distinct=true Frobenius_line_count11=true time=%o\n",TargetDegree,Cputime(tt);
load "P67/data/quotient45/target_places6_jbase_data.m";
assert #TargetPlaceKeys eq 11;kk:=[* *];ss:=[* *];for key in TargetPlaceKeys do z,sz:=ValidateKey(key);Append(~kk,z);Append(~ss,sz);end for;assert kk eq TargetPlaceKeys and #Seqset([k:k in kk])eq 11;tp:=[Index(kk,k):k in ss];assert tp eq TargetSigmaPermutation and Sort(tp)eq[1..11];assert &and[key[1]eq TargetJPolynomials[6]:key in kk];Append(~finitecounts,<TargetDegree,2^TargetDegree mod 67,11>);PrintFile(out,Sprintf("Append(~ValidatedTargetKeys,%o);Append(~ValidatedTargetPermutations,%o);",kk,tp));printf "NONIC67_FRESH_TARGET_KEYS_COMPLETE fiber=6 degree=%o count11=true distinct=true Frobenius_line_count11=true time=%o\n",TargetDegree,Cputime(tt);
load "P67/data/quotient45/target_places7_jbase_data.m";
assert #TargetPlaceKeys eq 11;kk:=[* *];ss:=[* *];for key in TargetPlaceKeys do z,sz:=ValidateKey(key);Append(~kk,z);Append(~ss,sz);end for;assert kk eq TargetPlaceKeys and #Seqset([k:k in kk])eq 11;tp:=[Index(kk,k):k in ss];assert tp eq TargetSigmaPermutation and Sort(tp)eq[1..11];assert &and[key[1]eq TargetJPolynomials[7]:key in kk];Append(~finitecounts,<TargetDegree,2^TargetDegree mod 67,11>);PrintFile(out,Sprintf("Append(~ValidatedTargetKeys,%o);Append(~ValidatedTargetPermutations,%o);",kk,tp));printf "NONIC67_FRESH_TARGET_KEYS_COMPLETE fiber=7 degree=%o count11=true distinct=true Frobenius_line_count11=true time=%o\n",TargetDegree,Cputime(tt);
PrintFile(out,Sprintf("FrobeniusFiberCounts:=%o;Complete:=true;",finitecounts));print "NONIC67_STANDALONE_FINITE_DATA_PRODUCER_COMPLETE sources55=true targets77=true all_fibers_independently_complete=true";quit;
