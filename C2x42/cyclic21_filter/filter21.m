
/* Fresh degree9 four-map filter on all first-stage surviving support types. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(7500*10^6);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
Attach("common/class_rational_reduction.m");
load "common/class_group_direct_RR_checked_v4_fastfb.m";
load "common/native_group_generators.m";
load "common/pair_af_place_linear.m";
Y:=MDX1(21,GF(5));C:=Curve(Y);print "C242_D9_NATIVE_CYCLIC21_GENUS",Genus(Y);assert Genus(Y) eq 5;cs:=Cusps(Y);os:=CuspOrbitsQ(Y);
c0:=Support([o:o in os|Degree(o) eq 1][1])[1];print "C242_D9_CYCLIC21_SETUP_READY",Cputime(t0);
AF36,am36:=AlgorithmicFunctionField(FunctionField(C));base36:=FunctionFieldDivisor(Divisor(c0));coords36:=[am36(FunctionField(C).i):i in [1,2]];invs:=ClassGroupAbelianInvariants(AF36:ReductionDivisor:=base36,Proof:=true);FF36:=RationalExtensionRepresentation(AF36);MDProvenClassGroupGenerators(FF36);
G,_,phi36,witness36:=MDClassGroupDirectRRCheckedV4FastFB(AF36,"C2x42/logs/cyclic21_filter/RR.log");
m2:=map<DivisorGroup(C)->G|D:->phi36(FunctionFieldDivisor(D))>;

hg:=[m2(o-Degree(o)*Divisor(c0)):o in os];H:=sub<G|hg>;assert Invariants(H) eq [364];QH,qh:=quo<G|H>;
cc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];comps:=[*{QH!0},{QH|},{QH|},{QH|},{QH|},{QH|}*];
for i in [1..#cs] do w:=Degree(cs[i]);if w le 5 then for d in [w..5] do comps[d+1] join:={z+qh(cc[i]):z in comps[d-w+1]};end for;end if;end for;
groupdata:=[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(v):v in hg],[Degree(c):c in cs],[Eltseq(v):v in cc],[[Eltseq(v):v in co]:co in comps]*];
print "C242_D9_CYCLIC21_GROUP_READY",Invariants(G),Invariants(H),[#s:s in comps],Cputime(t0);
Write("C2x42/data/cyclic21_filter/setup_group_arithmetic.txt",Sprint(groupdata):Overwrite:=true);
print "C242_D9_CYCLIC21_READY",Cputime(t0);

fp21:=FrobeniusPolynomial(JOne(21),5);assert &*[v:v in Invariants(G)|v ne 0] eq Integers()!Evaluate(fp21,1) and Integers()!Evaluate(fp21,1) eq 2184;
Write("C2x42/data/cyclic21_filter/quotient_coordinates.txt",Sprint([*Invariants(QH),[Eltseq(r):r in Rows(RelationMatrix(QH))],[Eltseq(qh(G.i)):i in [1..Ngens(G)]]*]):Overwrite:=true);
input:=eval Read("C2x42/data/filter_quotient14_linear_arithmetic.txt");singles:=input[3];small:=input[4];fields:=input[2];survivors:=(eval Read("C2x42/data/genus5_quotient/filter_919/filter_arithmetic.txt"))[10];assert #survivors eq 119;
own:=AssociativeArray();sm:=AssociativeArray();for d in [4..9] do own[d]:=[r:r in singles|r[1] eq d];end for;for d in [4,5] do sm[d]:=[r:r in small|r[1] eq d];end for;
R0:=PolynomialRing(GF(5));fieldobjects:=AssociativeArray();for a in fields do fieldobjects[a[1]]:=ext<GF(5)|R0!a[2]>;end for;
cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
procedure Image21(~ans,~cache,~images,row,mapno)
 d:=row[1];ky:=row[3];cachekey:=[mapno,d] cat ky;
 if IsDefined(cache,cachekey) then ans:=cache[cachekey];return;end if;
 K:=fieldobjects[d];z:=K.1;b:=Evaluate(R0!ky[1..d],z);c:=Evaluate(R0!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Order(Q) eq 42 and P ne 21*Q;Q21:=2*Q;
 if mapno eq 1 then EG:=E;QG:=Q21;
 else T:=mapno eq 2 select P else (mapno eq 3 select 21*Q else P+21*Q);EG,iso:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);QG:=iso(Q21);end if;
 assert Order(QG) eq 21;y:=C48PairAFPlaceLinear(coords36,TateNormalForm_xy(EG,QG));e:=Degree(y);assert d mod e eq 0;w:=d div e;ans:=w*phi36(Divisor(y)-e*base36);cache[cachekey]:=ans;
 Append(~images,[*cachekey,e,w,Eltseq(ans)*]);
end procedure;
for j in [1..#survivors] do
 desc:=survivors[j];sp:=[*own[desc[2]][desc[3]]*];if desc[1] ge 2 then Append(~sp,sm[desc[4]][desc[5]]);end if;ncdeg:=&+[r[1]:r in sp];cd:=9-ncdeg;assert cd in [0..5];checks:=[* *];passed:=true;
 for mapno in [1..4] do
  cl:=G!0;for row in sp do a:=G!0;Image21(~a,~cache,~images,row,mapno);cl+:=a;end for;
  yes:=-qh(cl) in comps[cd+1];Append(~checks,[*mapno,Eltseq(cl),yes*]);if not yes then passed:=false;break;end if;
 end for;Append(~tests,[*j,desc,cd,checks,passed*]);if passed then Append(~remaining,desc);end if;
 if j mod 20 eq 0 then print "C242_CYCLIC21_FILTER_PROGRESS",j,#remaining,Cputime(t0);end if;
end for;
Write("C2x42/data/cyclic21_filter/filter_arithmetic.txt",Sprint([*groupdata,images,tests,remaining*]):Overwrite:=true);fprintf witness36,"]\n";Flush(witness36);
print "C242_CYCLIC21_FILTER119_DONE",#tests,#remaining,"complete=true",Cputime(t0);quit;
