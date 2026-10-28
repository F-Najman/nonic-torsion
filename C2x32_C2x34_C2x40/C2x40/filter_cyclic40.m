/* Fresh degree9 four-map filter on all first-stage surviving support types. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
Attach("common/class_rational_reduction.m");
load "common/class_group_direct_RR_checked_v4_fastfb.m";
load "common/native_group_generators.m";
load "common/pair_af_place_linear.m";
Y:=MDX1(40,GF(3));C:=Curve(Y);print "C240_D9_NATIVE_CYCLIC40_GENUS",Genus(Y);assert Genus(Y) eq 25;cs:=Cusps(Y);os:=CuspOrbitsQ(Y);
c0:=Support([o:o in os|Degree(o) eq 1][1])[1];print "C240_D9_CYCLIC40_SETUP_READY",Cputime(t0);
AF36,am36:=AlgorithmicFunctionField(FunctionField(C));base36:=FunctionFieldDivisor(Divisor(c0));coords36:=[am36(FunctionField(C).i):i in [1,2]];invs:=ClassGroupAbelianInvariants(AF36:ReductionDivisor:=base36,Proof:=true);FF36:=RationalExtensionRepresentation(AF36);MDProvenClassGroupGenerators(FF36);
G,_,phi36,witness36:=MDClassGroupDirectRRCheckedV4FastFB(AF36,"C2x32_C2x34_C2x40/logs/C2x40/filter_cyclic40_RR_identity.log");
m2:=map<DivisorGroup(C)->G|D:->phi36(FunctionFieldDivisor(D))>;

hg:=[m2(o-Degree(o)*Divisor(c0)):o in os];H:=sub<G|hg>;assert Invariants(H) eq [2,2,2,8,120,895440];QH,qh:=quo<G|H>;
cc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];comps:=[*{QH!0},{QH|},{QH|},{QH|}*];
for i in [1..#cs] do w:=Degree(cs[i]);if w le 3 then for d in [w..3] do comps[d+1] join:={z+qh(cc[i]):z in comps[d-w+1]};end for;end if;end for;
groupdata:=[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(v):v in hg],[Degree(c):c in cs],[Eltseq(v):v in cc],[[Eltseq(v):v in co]:co in comps]*];
print "C240_D9_CYCLIC40_GROUP_READY",Invariants(G),Invariants(H),[#s:s in comps],Cputime(t0);
data:=eval Read("C2x32_C2x34_C2x40/data/C2x40/filter_quotient20_linear_arithmetic.txt");singles:=data[3];small:=data[4];survivors:=data[7];fields:=data[2];assert #singles eq 106 and #data[5] eq 128;
survivors:=(eval Read("C2x32_C2x34_C2x40/data/C2x40/filter_quotient20_isogeny_arithmetic.txt"))[4];assert #survivors eq 27;
own:=AssociativeArray();sm:=AssociativeArray();for d in [4..9] do own[d]:=[r:r in singles|r[1] eq d];end for;for d in [4,5] do sm[d]:=[r:r in small|r[1] eq d];end for;
R:=PolynomialRing(GF(3));fieldobjects:=AssociativeArray();for a in fields do fieldobjects[a[1]]:=ext<GF(3)|R!a[2]>;end for;
cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
procedure Push(~answer,~cache,~images,row,mapno,Y,c0,m2)
 D:=row[1];ky:=row[3];cachekey:=[mapno,D] cat ky;if IsDefined(cache,cachekey) then answer:=cache[cachekey];return;end if;
 K:=fieldobjects[D];z:=K.1;b:=Evaluate(R!ky[1..D],z);c:=Evaluate(R!ky[D+1..2*D],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];
 assert Order(Q) eq 40 and P ne 20*Q;
 if mapno le 2 then EG:=E;QG:=mapno eq 1 select Q else P+Q;
 else T:=mapno eq 3 select P else P+20*Q;EG,phi:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);QG:=phi(Q);end if;assert Order(QG) eq 40;
 y:=C48PairAFPlaceLinear(coords36,TateNormalForm_xy(EG,QG));e:=Degree(y);assert D mod e eq 0;w:=D div e;
 answer:=w*phi36(Divisor(y)-e*base36);cache[cachekey]:=answer;Append(~images,[*cachekey,e,w,Eltseq(answer)*]);
end procedure;
for k in [1..#survivors] do
 desc:=survivors[k];sp:=[*own[desc[2]][desc[3]]*];
 if desc[1] ge 2 then Append(~sp,sm[desc[4]][desc[5]]);end if;
 if desc[1] eq 3 then Append(~sp,sm[desc[6]][desc[7]]);end if;
 ncdeg:=&+[r[1]:r in sp];cd:=9-ncdeg;assert cd in [0..3];passed:=true;checks:=[* *];
 for mapno in [1..4] do
  cl:=G!0;for r in sp do a:=G!0;Push(~a,~cache,~images,r,mapno,Y,c0,m2);cl+:=a;end for;
  yes:=-qh(cl) in comps[cd+1];Append(~checks,[*mapno,Eltseq(cl),yes*]);if not yes then passed:=false;break;end if;
 end for;
 Append(~tests,[*k,desc,cd,checks,passed*]);if passed then Append(~remaining,desc);end if;
 if k mod 1 eq 0 then
  Write("C2x32_C2x34_C2x40/data/C2x40/filter_cyclic40_arithmetic.txt",Sprint([*groupdata,images,tests,remaining*]):Overwrite:=true);
  print "C240_D9_CYCLIC40_PROGRESS",k,#remaining,Cputime(t0);
 end if;
end for;
Write("C2x32_C2x34_C2x40/data/C2x40/filter_cyclic40_arithmetic.txt",Sprint([*groupdata,images,tests,remaining*]):Overwrite:=true);fprintf witness36,"]\n";Flush(witness36);
printf "C240_D9_CYCLIC40_DONE tested=%o remaining=%o complete=true cpu=%o\n",#tests,#remaining,Cputime(t0);quit;
