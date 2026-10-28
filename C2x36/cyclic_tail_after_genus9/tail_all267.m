/* Fresh degree9 four-map filter on all first-stage surviving support types.
   Rerun made for the repository (27 October 2026): tail.m applied to all 267 genus-9 survivors
   exported by the diagonal quotient filter, in one fail-closed batch run (tail.m itself treated
   the 200 survivors not covered by the interrupted run of filter_cyclic36_linear.m). Changes:
   the list of survivors, its size, and the output paths. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
Attach("common/class_rational_reduction.m");
load "common/class_group_direct_RR_checked_v4_fastfb.m";
load "common/native_group_generators.m";
load "common/pair_af_place_linear.m";
Y:=MDX1(36,GF(5));C:=Curve(Y);cs:=Cusps(Y);os:=CuspOrbitsQ(Y);assert Genus(Y) eq 17;
c0:=Support([o:o in os|Degree(o) eq 1][1])[1];print "C236_D9_CYCLIC36_SETUP_READY",Cputime(t0);
AF36,am36:=AlgorithmicFunctionField(FunctionField(C));base36:=FunctionFieldDivisor(Divisor(c0));coords36:=[am36(FunctionField(C).i):i in [1,2]];invs:=ClassGroupAbelianInvariants(AF36:ReductionDivisor:=base36,Proof:=true);FF36:=RationalExtensionRepresentation(AF36);MDProvenClassGroupGenerators(FF36);
G,_,phi36,witness36:=MDClassGroupDirectRRCheckedV4FastFB(AF36,"C2x36/logs/cyclic_tail_after_genus9/tail_all267_RR_identity.log");
m2:=map<DivisorGroup(C)->G|D:->phi36(FunctionFieldDivisor(D))>;
assert Invariants(G) eq [4,4,12,12,36,7812,7812,0];
hg:=[m2(o-Degree(o)*Divisor(c0)):o in os];H:=sub<G|hg>;assert Invariants(H) eq [12,252,7812];QH,qh:=quo<G|H>;assert &*[v:v in Invariants(QH)|v ne 0] eq 214272;
cc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];comps:=[*{QH!0},{QH|},{QH|},{QH|}*];
for i in [1..#cs] do w:=Degree(cs[i]);if w le 3 then for d in [w..3] do comps[d+1] join:={z+qh(cc[i]):z in comps[d-w+1]};end for;end if;end for;
groupdata:=[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(v):v in hg],[Degree(c):c in cs],[Eltseq(v):v in cc],[[Eltseq(v):v in co]:co in comps]*];
print "C236_D9_CYCLIC36_GROUP_READY",Invariants(G),Invariants(H),[#s:s in comps],Cputime(t0);
data:=eval Read("C2x36/data/filter_quotient18_linear_arithmetic.txt");singles:=data[3];small:=data[4];survivors:=data[7];fields:=data[2];assert #singles eq 11800 and #data[5] eq 11784 and #data[6] eq 2352;
survivors:=(eval Read("C2x36/data/diagonal_quotient/quotient_filter_arithmetic.txt"))[4];assert #survivors eq 267;
own:=AssociativeArray();sm:=AssociativeArray();for d in [3..9] do own[d]:=[r:r in singles|r[1] eq d];end for;for d in [3..6] do sm[d]:=[r:r in small|r[1] eq d];end for;
R:=PolynomialRing(GF(5));fieldobjects:=AssociativeArray();for a in fields do fieldobjects[a[1]]:=ext<GF(5)|R!a[2]>;end for;
cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
procedure Push(~answer,~cache,~images,row,mapno,Y,c0,m2)
 D:=row[1];ky:=row[3];cachekey:=[mapno,D] cat ky;if IsDefined(cache,cachekey) then answer:=cache[cachekey];return;end if;
 K:=fieldobjects[D];z:=K.1;b:=Evaluate(R!ky[1..D],z);c:=Evaluate(R!ky[D+1..2*D],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];
 assert Order(Q) eq 36 and P ne 18*Q;
 if mapno le 2 then EG:=E;QG:=mapno eq 1 select Q else P+Q;
 else T:=mapno eq 3 select P else P+18*Q;EG,phi:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);QG:=phi(Q);end if;assert Order(QG) eq 36;
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
 if k mod 20 eq 0 then
  Write("C2x36/data/cyclic_tail_after_genus9/tail_all267_arithmetic.txt",Sprint([*groupdata,images,tests,remaining*]):Overwrite:=true);
  print "C236_D9_CYCLIC36_PROGRESS",k,#remaining,Cputime(t0);
 end if;
end for;
Write("C2x36/data/cyclic_tail_after_genus9/tail_all267_arithmetic.txt",Sprint([*groupdata,images,tests,remaining*]):Overwrite:=true);fprintf witness36,"]\n";Flush(witness36);
printf "C236_D9_REDUCED_CYCLIC36_TAIL_DONE tested=%o remaining=%o complete=true cpu=%o\n",#tests,#remaining,Cputime(t0);quit;
