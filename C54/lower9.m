SetColumns(0);SetQuitOnError(true);tt:=Cputime();
load "C54/preflight_module9.m";
load "common/pair_af_place_linear.m";
FF5:=FunctionField(C5);AF5,am5:=AlgorithmicFunctionField(FF5);AF5:=RationalExtensionRepresentation(AF5);
GAF,_,phi:=ClassGroup(AF5:Proof:=true);assert GAF cmpeq G5;
coords:=[AF5!am5(FF5.i):i in [1..2]];A0:=FunctionFieldDivisor(Divisor(c05));assert Degree(A0) eq 1;
assert &and[&and[a in H5:a in comps[k]]:k in [0..6]];
function NativePlaceClass(P)
 values:=[Evaluate(FF5.i,P):i in [1..2]];Q:=C48PairAFPlaceLinear(coords,values);assert Degree(Q) eq Degree(P);
 return phi(Divisor(Q)-Degree(Q)*A0);
end function;
function Representatives(ps)
 seen:=Seqset(ps);answer:=[];
 while #seen gt 0 do
  P:=Rep(seen);orbit:=Seqset(DiamondOrbit(X5,P));assert orbit subset seen;
  Append(~answer,Index(ps,P));seen diff:=orbit;
 end while;return answer;
end function;
allps:=[[]:d in [1..6]];good:=[[]:d in [1..6]];inert:=[[]:d in [1..6]];r2counts:=[];
for d in [2..6] do
 allps[d]:=NoncuspidalPlaces(X5,d);nrlist:=[];
 for P in allps[d] do
  E:=EllipticCurve(X5,P);b2,b4,b6,b8:=Explode(bInvariants(E));RR<T>:=PolynomialRing(BaseRing(E));nr:=#Roots(4*T^3+b2*T^2+2*b4*T+b6);assert nr in [0,1,3];Append(~nrlist,nr);
  if nr gt 0 then Append(~good[d],P);else Append(~inert[d],P);end if;
 end for;
 Append(~r2counts,nrlist);printf "NONIC54_LOWER_POINTS degree=%o all=%o rational2=%o inert=%o time=%o\n",d,#allps[d],#good[d],#inert[d],Cputime(tt);
end for;
assert [#allps[d]:d in [2..6]] eq [9,57,171,612,2700] and [#good[d]:d in [2..6]] eq [0,30,171,396,1890];
classes:=AssociativeArray();reps:=AssociativeArray();
for d in [3..6] do
 reps[d]:=Representatives(good[d]);classes[d]:=[];
 for P in good[d] do Append(~classes[d],NativePlaceClass(P));end for;
 printf "NONIC54_LOWER_CLASSES degree=%o places=%o diamond_reps=%o time=%o\n",d,#good[d],#reps[d],Cputime(tt);
end for;
out:=Open("C54/data/lower_testclasses.py","w");
fprintf out,"relations2=%o\nrelations5=%o\ncusp2=%o\ncusp5=%o\ncandidate2=%o\ntests=[\n",[Eltseq(r):r in Rows(RelationMatrix(G2))],[Eltseq(r):r in Rows(RelationMatrix(G5))],[Eltseq(v):v in qc2],[Eltseq(v):v in qc5],[Eltseq(v):v in candidate2];Flush(out);
checked:=0;members:=0;passing:=0;
procedure Test(a,d,desc,~checked,~members,~passing)
 for cc in comps[9-d] do
  aa:=a+cc;member:=aa in H5;passes:=false;im2:=[];preimage:=[];
  if member then
   members+:=1;v:=aa@@h5;assert h5(v) eq aa;im2:=Eltseq(h2(v));preimage:=Eltseq(v);if h2(v) in candidate2 then passes:=true;passing+:=1;end if;
  end if;
  checked+:=1;fprintf out,"[%o,%o,%o,%o,%o,%o],\n",desc,Eltseq(aa),member select 1 else 0,passes select 1 else 0,im2,preimage;
 end for;
end procedure;
for d in [3..6] do
 for i in reps[d] do Test(classes[d][i],d,[d,i],~checked,~members,~passing);end for;
 printf "NONIC54_LOWER_SINGLE_DONE degree=%o checked=%o members=%o passing=%o time=%o\n",d,checked,members,passing,Cputime(tt);
end for;
for d in [2,3] do for i in Representatives(inert[d]) do Test(3*NativePlaceClass(inert[d][i]),3*d,[d,d,d,i],~checked,~members,~passing);end for;end for;
for de in [<3,3>,<3,4>,<3,5>,<3,6>,<4,4>,<4,5>] do
 d,e:=Explode(de);
 for i in reps[d],j in [1..#classes[e]] do Test(classes[d][i]+classes[e][j],d+e,[d,e,i,j],~checked,~members,~passing);end for;
 printf "NONIC54_LOWER_MIXED_DONE degrees=%o checked=%o members=%o passing=%o time=%o\n",de,checked,members,passing,Cputime(tt);
end for;
for i in reps[3],j in [1..#classes[3]],k in [j..#classes[3]] do Test(classes[3][i]+classes[3][j]+classes[3][k],9,[3,3,3,i,j,k],~checked,~members,~passing);end for;
printf "NONIC54_LOWER_TRIPLES_DONE checked=%o members=%o passing=%o time=%o\n",checked,members,passing,Cputime(tt);
data7:=eval Read("C54/data/degree8/good7_pointdata.txt");def7,records7:=Explode(data7);assert #records7 eq 828;
R7:=PolynomialRing(GF(5));f7:=R7!def7;assert Degree(f7) eq 7 and IsIrreducible(f7);K7<a7>:=ext<GF(5)|f7>;weighted7:=0;
for i in [1..#records7] do
 xyco,n,r2,jco:=Explode(records7[i]);assert n mod 54 eq 0 and r2 in [1,3];xy:=[Evaluate(R7!c,a7):c in xyco];assert LCM([Degree(MinimalPolynomial(z,GF(5))):z in xy]) eq 7;
 Q:=C48PairAFPlaceLinear(coords,xy);assert Degree(Q) eq 7;aa:=phi(Divisor(Q)-7*A0);Test(aa,7,[7,i],~checked,~members,~passing);weighted7+:=9*r2;
 if i mod 100 eq 0 then printf "NONIC54_GOOD7_PROGRESS tested=%o checked=%o members=%o passing=%o time=%o\n",i,checked,members,passing,Cputime(tt);end if;
end for;
assert weighted7 eq 11412;
fprintf out,"]\ncomplete=True\n";delete out;
printf "NONIC54_LOWER_DONE good7=828 weighted7=11412 checked=%o members=%o passing=%o complete=true time=%o\n",checked,members,passing,Cputime(tt);quit;
