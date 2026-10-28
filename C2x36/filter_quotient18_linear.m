/* New exact degree9 necessary filter; no earlier-degree exclusion is used. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/class_rational_reduction.m");
load "common/class_group_direct_RR_checked_v4_fastfb.m";
load "common/pair_af_place_linear.m";
Y:=MDX11(2,18,GF(5));C:=Curve(Y);assert Genus(C) eq 7;cs:=Cusps(Y);c0:=Rep({c:c in cs|Degree(c) eq 1});
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
AF18,am18:=AlgorithmicFunctionField(FunctionField(C));coords18:=[am18(FunctionField(C)!f):f in Y`_coordinates];base18:=FunctionFieldDivisor(Divisor(c0));Gfast,_,phi18,witness18:=MDClassGroupDirectRRCheckedV4FastFB(AF18,"C2x36/logs/quotient18_linear_RR_identity.log");assert Gfast cmpeq G;
m2:=map<DivisorGroup(C)->G|D:->phi18(FunctionFieldDivisor(D))>;
cc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];hg:=cc;H:=sub<G|hg>;
/* Published rational cusp order equals this finite cusp order.
   Its reduction lies in H and injects, so H is exactly the rational image. */
assert Invariants(G) eq [6,84,252,0] and Invariants(H) eq [2,42,126];
groupdata:=[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Degree(c):c in cs],[Eltseq(v):v in cc],[Eltseq(v):v in hg],Invariants(H)*];
function BcKey(bc,D)
 keys:=[];for i in [1..D] do Append(~keys,&cat[[Integers()!x:x in Eltseq(z,GF(5))]:z in bc]);bc:=[z^5:z in bc];end for;return Sort(keys)[1];
end function;
function BC(E,P,Q)
 ai:=aInvariants(E);xx:=Q[1]-P[1];yy:=Q[2];return [xx^3/yy^2,(ai[2]+3*P[1])*xx^2/yy^2];
end function;
procedure PushClass(~cache,~answer,~edata,b,c,D,Y,c0,m2)
 dd:=b*(1-b-c);s:=(3*b^2+2*c*b+dd)/(2*b);xx:=s^2-c-2*b;yy:=-b+s*(b-xx);assert xx ne 0 and yy ne 0;
 bc:=[xx^3/yy^2,c*xx^2/yy^2];ky:=[D] cat BcKey(bc,D);
 if IsDefined(cache,ky) then answer:=cache[ky][1];edata:=cache[ky][2];return;end if;
 y:=C48PairAFPlaceLinear(coords18,bc);e:=Degree(y);assert D mod e eq 0 and D div e in [1,2,4];
 w:=D div e;answer:=w*phi18(Divisor(y)-e*base18);edata:=[e,w];cache[ky]:=<answer,edata>;
end procedure;
R:=PolynomialRing(GF(5));fields:=[* *];single:=[* *];small:=[* *];pairs:=[* *];triples:=[* *];survivors:=[* *];own:=AssociativeArray();cache:=AssociativeArray();smalls:=AssociativeArray();coverage:=[* *];
for D in [3..9] do
 a:=eval Read(Sprintf("C2x36/data/M2N36_D9_direct_full%o_arithmetic.txt",D));assert a[1] eq D and a[3];def:=a[2];records:=a[9];
 K<z>:=ext<GF(5)|R!def>;Append(~fields,[*D,def*]);own[D]:=[* *];passed:=0;allkeys:={};orbcov:=[* *];
 for inputi in [1..#records] do
  d,key:=Explode(records[inputi]);if d ne D then continue;end if;i:=#own[D]+1;
  b:=Evaluate(R!key[1..D],z);c:=Evaluate(R!key[D+1..2*D],z);cl:=G!0;ed:=[];
  PushClass(~cache,~cl,~ed,b,c,D,Y,c0,m2);yes:=cl in H;row:=[*D,i,key,ed[1],ed[2],Eltseq(cl),yes select 1 else 0*];Append(~single,row);Append(~own[D],row);
  if yes then passed+:=1;Append(~survivors,[*1,D,i*]);end if;
  if D le 6 then
   E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P:=E![0,0];Q:=E![b,b];assert Order(P) eq 2 and Order(Q) eq 36 and P ne 18*Q;
   keys:={BcKey(BC(E,P+18*bb*Q,cc*P+aa*Q),D):bb,cc in [0,1],aa in [1,5,7,11,13,17]};assert IsDisjoint(allkeys,keys);allkeys join:=keys;Append(~orbcov,[*i,Setseq(keys)*]);
  end if;
  if #own[D] mod 500 eq 0 then print "C236_D9_QUOTIENT_SINGLE_PROGRESS",D,#own[D],passed,Cputime(t0);end if;
 end for;
 print "C236_D9_QUOTIENT_SINGLE_TYPE",D,#own[D],passed,Cputime(t0);
 if D le 6 then
  smalls[D]:=[* *];keys:=Sort(Setseq(allkeys));assert #keys eq [48,120,720,2904][D-2];Append(~coverage,[*D,orbcov,keys*]);
  for j in [1..#keys] do
   ky:=keys[j];b:=Evaluate(R!ky[1..D],z);c:=Evaluate(R!ky[D+1..2*D],z);cl:=G!0;ed:=[];PushClass(~cache,~cl,~ed,b,c,D,Y,c0,m2);
   row:=[*D,j,ky,ed[1],ed[2],Eltseq(cl)*];Append(~small,row);Append(~smalls[D],row);
  end for;
  print "C236_D9_QUOTIENT_ALL_SMALL",D,#smalls[D],Cputime(t0);
 end if;
 Write("C2x36/data/filter_quotient18_linear_arithmetic.txt",Sprint([*groupdata,fields,single,small,pairs,triples,survivors,coverage*]):Overwrite:=true);
end for;
assert #single eq 11800;
for typ in [[3,3],[3,4],[3,5],[3,6],[4,4],[4,5]] do
 D,E:=Explode(typ);passed:=0;tested:=0;
 for i in [1..#own[D]] do
  cl1:=G!own[D][i][6];
  for j in [1..#smalls[E]] do
   cl:=cl1+G!smalls[E][j][6];yes:=cl in H;tested+:=1;Append(~pairs,[*D,i,E,j,Eltseq(cl),yes select 1 else 0*]);
   if yes then passed+:=1;Append(~survivors,[*2,D,i,E,j*]);end if;
  end for;
 end for;
 print "C236_D9_QUOTIENT_PAIR_TYPE",D,E,tested,passed,Cputime(t0);
end for;
for i in [1..#own[3]] do
 cl1:=G!own[3][i][6];
 for j in [1..#smalls[3]],k in [j..#smalls[3]] do
  cl:=cl1+G!smalls[3][j][6]+G!smalls[3][k][6];yes:=cl in H;Append(~triples,[*3,i,j,k,Eltseq(cl),yes select 1 else 0*]);
  if yes then Append(~survivors,[*3,3,i,3,j,3,k*]);end if;
 end for;
end for;
Write("C2x36/data/filter_quotient18_linear_arithmetic.txt",Sprint([*groupdata,fields,single,small,pairs,triples,survivors,coverage*]):Overwrite:=true);
fprintf witness18,"]\n";Flush(witness18);
printf "C236_D9_QUOTIENT18_LINEAR_DONE single=%o pair=%o triple=%o remaining=%o complete=true cpu=%o\n",#single,#pairs,#triples,#survivors,Cputime(t0);quit;
