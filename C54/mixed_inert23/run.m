SetColumns(0);SetQuitOnError(true);SetMemoryLimit(6*10^9);tp:=Cputime();
load "C54/data/projection2_v2/inputs.m";
load "C54/projection2_v2/preflight_ordinary_snapshot.m";
assert sig2 eq ExpectedSigs;
FF5:=FunctionField(C5);AF5,am5:=AlgorithmicFunctionField(FF5);AF5:=RationalExtensionRepresentation(AF5);coords:=[AF5!am5(FF5.i):i in [1..2]];
G5AF,afinv,phi:=ClassGroup(AF5:Proof:=true);assert G5AF cmpeq G5;A0:=FunctionFieldDivisor(Divisor(c05));
load "common/pair_af_place_linear.m";
function PairClass(xy,d)
 Q:=C48PairAFPlaceLinear(coords,xy);e:=Degree(Q);assert d mod e eq 0;
 return phi((d div e)*Divisor(Q)-d*A0),e;
end function;
function ClassCondition(a)
 if a notin H5 then return false,[],[];end if;
 v:=a@@h5;assert h5(v) eq a;im:=h2(v);return im in candidate2,Eltseq(im),Eltseq(v);
end function;
function SecondClasses(E,P,d)
 assert 27*P eq E!0 and 9*P ne E!0;
 K:=BaseRing(E);R<T>:=PolynomialRing(K);bb2,bb4,bb6,bb8:=Explode(bInvariants(E));rr:=Roots(4*T^3+bb2*T^2+2*bb4*T+bb6);assert #rr in [1,3];
 out:=[* *];
 for ri in [1..#rr] do
  root:=rr[ri][1];ai:=aInvariants(E);T2:=E![root,-(ai[1]*root+ai[3])/2];assert T2 ne E!0 and 2*T2 eq E!0;
  EQ,iso:=IsogenyFromKernel(E,T-root);PQ:=iso(P);assert 27*PQ eq EQ!0 and 9*PQ ne EQ!0 and Degree(iso) eq 2;
  xy:=TateNormalForm_xy(EQ,PQ);aa,e:=PairClass(xy,d);Append(~out,<aa,e,ri>);
 end for;return out;
end function;

function Representatives(ps)
 seen:=Seqset(ps);answer:=[];
 while #seen gt 0 do P:=Rep(seen);orbit:=Seqset(DiamondOrbit(X5,P));assert orbit subset seen;Append(~answer,Index(ps,P));seen diff:=orbit;end while;return answer;
end function;
ps2:=NoncuspidalPlaces(X5,2);assert #ps2 eq 9;
good3:=[];for Q in NoncuspidalPlaces(X5,3) do E:=EllipticCurve(X5,Q);R:=PolynomialRing(BaseRing(E));b2e,b4e,b6e,b8e:=Explode(bInvariants(E));if #Roots(4*R.1^3+b2e*R.1^2+2*b4e*R.1+b6e) gt 0 then Append(~good3,Q);end if;end for;assert #good3 eq 30;
reps2:=Representatives(ps2);assert #reps2 eq 1;
out:=Open("C54/data/mixed_inert23/results.py","w");fprintf out,"relations2=%o\nrelations5=%o\ncusp2=%o\ncusp5=%o\ncandidate2=%o\ntests=[\n",[Eltseq(r):r in Rows(RelationMatrix(G2))],[Eltseq(r):r in Rows(RelationMatrix(G5))],[Eltseq(v):v in qc2],[Eltseq(v):v in qc5],[Eltseq(v):v in candidate2];Flush(out);
secondOut:=Open("C54/data/mixed_inert23/second_results.py","w");fprintf secondOut,"rows=[\n";tested:=0;members:=0;passes:=0;secondPasses:=0;
for i in reps2 do
 Q2:=ps2[i];E2:=EllipticCurve(X5,Q2);P2:=LevelStructure(X5,Q2)`P;R2<T>:=PolynomialRing(BaseRing(E2));bb2,bb4,bb6,bb8:=Explode(bInvariants(E2));poly:=4*T^3+bb2*T^2+2*bb4*T+bb6;assert IsIrreducible(poly) and #Roots(poly) eq 0;
 aa2,e2:=PairClass([Evaluate(FF5.k,Q2):k in [1..2]],2);assert e2 eq 2;
 for j in [1..#good3] do
  Q3:=good3[j];aa3,e3:=PairClass([Evaluate(FF5.k,Q3):k in [1..2]],3);assert e3 eq 3;aa:=3*aa2+aa3;pass,im,pre:=ClassCondition(aa);member:=aa in H5;if member then members+:=1;end if;if pass then passes+:=1;end if;tested+:=1;
  fprintf out,"[%o,%o,%o,%o,%o,%o],\n",[2,2,2,3,i,j],Eltseq(aa),member select 1 else 0,pass select 1 else 0,im,pre;Flush(out);
  if pass then
   L<root>:=ext<BaseRing(E2)|poly>;EL:=BaseChange(E2,L);PL:=EL![L!c:c in Eltseq(P2)];ss2:=SecondClasses(EL,PL,6);assert #ss2 eq 3 and #Seqset([v[1]:v in ss2]) eq 1;
   E3:=EllipticCurve(X5,Q3);P3:=LevelStructure(X5,Q3)`P;ss3:=SecondClasses(E3,P3,3);
   for v in ss3 do
    bb:=ss2[1][1]+v[1];pp,iim,ppr:=ClassCondition(bb);if pp then secondPasses+:=1;end if;
    fprintf secondOut,"[%o,%o,%o,%o,%o,%o],\n",[i,j,v[3]],Eltseq(bb),bb in H5 select 1 else 0,pp select 1 else 0,iim,ppr;Flush(secondOut);
   end for;
  end if;
 end for;
end for;
fprintf out,"]\ncomplete=True\n";delete out;fprintf secondOut,"]\ncomplete=True\n";delete secondOut;
printf "NONIC54_MIXED_INERT23_DONE representatives=%o good3=30 tests=%o members=%o first_passes=%o second_passes=%o complete=true time=%o\n",#reps2,tested,members,passes,secondPasses,Cputime(tp);quit;
