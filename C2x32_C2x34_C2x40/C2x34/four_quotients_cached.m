/* Complete fresh degree9 necessary four-map support filter. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");Attach("common/fast_diamond_hash_m2.m");Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("C2x32_C2x34_C2x40/C2x34/deps/fast_full_aut_hash_m2_34.m");
load "common/native_group_generators.m";
load "common/class_group_direct_RR_checked_v4_fastfb.m";
Y:=MDX1(34,GF(3));C:=Curve(Y);assert Genus(C) eq 21;os:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in CuspOrbitsQ(Y)];cs:=Cusps(Y);c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
AF:=AlgorithmicFunctionField(FunctionField(C));base:=FunctionFieldDivisor(Divisor(c0));invs:=ClassGroupAbelianInvariants(AF:ReductionDivisor:=base,Proof:=true);F:=RationalExtensionRepresentation(AF);MDProvenClassGroupGenerators(F);
G,_,phi,witness:=MDClassGroupDirectRRCheckedV4FastFB(AF,"C2x32_C2x34_C2x40/logs/C2x34/class_cached_RR_identity.log");
hg:=[phi(FunctionFieldDivisor(o-Degree(o)*Divisor(c0))):o in os];H:=sub<G|hg>;assert Invariants(H) eq [8760,595680];GG,pi:=quo<G|H>;assert &*[v:v in Invariants(GG)|v ne 0] eq 128;
cc:=[phi(FunctionFieldDivisor(Divisor(c)-Degree(c)*Divisor(c0))):c in cs];allowed:=AssociativeArray();for d in [0..5] do allowed[d]:={GG|};end for;allowed[0]:={GG!0};
for i in [1..#cs] do w:=Degree(cs[i]);if w le 5 then for d in [w..5] do allowed[d] join:={z+pi(cc[i]):z in allowed[d-w]};end for;end if;end for;
groupdata:=[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(g):g in hg],[Degree(c):c in cs],[Eltseq(a):a in cc]*];
print "C234_D9_FOURQUOTIENT_GROUP",Invariants(G),Invariants(H),[#allowed[d]:d in [0..5]],Cputime(tt);
src:=eval Read("C2x32_C2x34_C2x40/data/C2x34/source_coverage.txt");source,smalldata,coverage:=Explode(src);
fr:=eval Read("C2x32_C2x34_C2x40/data/C2x34/source_frobenius.txt");fp:=PolynomialRing(Integers())!fr[1];xcsdegrees:=fr[2];nc:=fr[3];assert nc eq [0,0,0,36,48,48,144,738,2160];
R:=PolynomialRing(GF(3));own:=AssociativeArray();small:=AssociativeArray();fields:=AssociativeArray();
for d in [4..9] do own[d]:=[r:r in source|r[1] eq d];co:=[r:r in coverage|r[1] eq d][1];fields[d]:=ext<GF(3)|R!co[2]>;if d le 5 then small[d]:=[r:r in smalldata|r[1] eq d];end if;end for;
cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
procedure Push(~answer,~cache,~images,row,variant)
 d:=row[1];ky:=row[4];cachekey:=[variant,d] cat ky;if IsDefined(cache,cachekey) then answer:=cache[cachekey];return;end if;
 K:=fields[d];z:=K.1;b:=Evaluate(R!ky[1..d],z);c:=Evaluate(R!ky[d+1..2*d],z);E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P:=E![0,0];Q:=E![b,b];assert Order(Q) eq 34 and P ne 17*Q;
 if variant le 2 then EG:=E;U:=variant eq 1 select Q else P+Q;else T:=variant eq 3 select P else P+17*Q;EG,ph:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);U:=ph(Q);end if;assert Order(U) eq 34;
 pair:=TateNormalForm_xy(EG,U);FC:=FunctionField(C);yp:=MDPlaceCommonZeros([FC.1,FC.2],pair);e:=Degree(yp);assert d mod e eq 0;w:=d div e;answer:=w*phi(FunctionFieldDivisor(Divisor(yp)-e*Divisor(c0)));cache[cachekey]:=answer;Append(~images,[*cachekey,e,w,Eltseq(answer)*]);
end procedure;
supports:=[*[1,r[1],r[2]]:r in source*];for ty in [[4,4],[4,5]] do for i in [1..#own[ty[1]]],j in [1..#small[ty[2]]] do Append(~supports,[2,ty[1],i,ty[2],j]);end for;end for;
for k in [1..#supports] do
 desc:=supports[k];sp:=[*own[desc[2]][desc[3]]*];if desc[1] eq 2 then Append(~sp,small[desc[4]][desc[5]]);end if;cd:=9-&+[r[1]:r in sp];assert cd in [0..5];passed:=true;checks:=[* *];
 for variant in [1..4] do cl:=G!0;for row in sp do a:=G!0;Push(~a,~cache,~images,row,variant);cl+:=a;end for;flag:=-pi(cl) in allowed[cd];Append(~checks,[*variant,Eltseq(cl),flag*]);if not flag then passed:=false;break;end if;end for;
 Append(~tests,[*k,desc,cd,checks,passed*]);if passed then Append(~remaining,desc);end if;
 if k mod 20 eq 0 then Write("C2x32_C2x34_C2x40/data/C2x34/four_quotients_cached_arithmetic.txt",Sprint([*groupdata,source,smalldata,coverage,images,tests,remaining,Coefficients(fp),xcsdegrees*]):Overwrite:=true);print "C234_D9_FOURQUOTIENT_PROGRESS",k,#remaining,Cputime(tt);end if;
end for;
Write("C2x32_C2x34_C2x40/data/C2x34/four_quotients_cached_arithmetic.txt",Sprint([*groupdata,source,smalldata,coverage,images,tests,remaining,Coefficients(fp),xcsdegrees*]):Overwrite:=true);fprintf witness,"]\n";Flush(witness);
printf "C234_D9_FOURQUOTIENT_CACHED_DONE tested=%o remaining=%o complete=true cpu=%o\n",#tests,#remaining,Cputime(tt);quit;
