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
out:=Open("C54/data/projection2_v2/results.py","w");
fprintf out,"relations2=%o\nrelations5=%o\ncusp2=%o\ncusp5=%o\ncandidate2=%o\nrows=[\n",[Eltseq(r):r in Rows(RelationMatrix(G2))],[Eltseq(r):r in Rows(RelationMatrix(G5))],[Eltseq(v):v in qc2],[Eltseq(v):v in qc5],[Eltseq(v):v in candidate2];Flush(out);
good:=AssociativeArray();for d in [3..6] do
 good[d]:=[];for Q in NoncuspidalPlaces(X5,d) do E:=EllipticCurve(X5,Q);b2e,b4e,b6e,b8e:=Explode(bInvariants(E));rr:=PolynomialRing(BaseRing(E));if #Roots(4*rr.1^3+b2e*rr.1^2+2*b4e*rr.1+b6e) gt 0 then Append(~good[d],Q);end if;end for;
end for;
assert [#good[d]:d in [3..6]] eq [30,171,396,1890];
for k in [1..#LowerDescriptors] do
 d,e,i,j:=Explode(LowerDescriptors[k]);Q1:=good[d][i];Q2:=good[e][j];E1:=EllipticCurve(X5,Q1);P1:=LevelStructure(X5,Q1)`P;E2:=EllipticCurve(X5,Q2);P2:=LevelStructure(X5,Q2)`P;
 a1,ee1:=PairClass([Evaluate(FF5.i,Q1):i in [1..2]],d);a2,ee2:=PairClass([Evaluate(FF5.i,Q2):i in [1..2]],e);assert a1+a2 eq h5(L!LowerPreimages[k]);
 ss1:=SecondClasses(E1,P1,d);ss2:=SecondClasses(E2,P2,e);
 for s1 in ss1,s2 in ss2 do
  aa:=s1[1]+s2[1];pass,im,pre:=ClassCondition(aa);
  fprintf out,"[%o,%o,%o,%o,%o,%o,%o,%o],\n",[0,k,s1[3],s2[3]],[d,e],[s1[2],s2[2]],Eltseq(aa),aa in H5 select 1 else 0,pass select 1 else 0,im,pre;Flush(out);
 end for;
 print "NONIC54_SECOND_LOWER",k,#ss1,#ss2,Cputime(tp);
end for;
data:=eval Read("C54/data/good9_pointdata.txt");def,records:=Explode(data);Rp:=PolynomialRing(GF(5));K9<a9>:=ext<GF(5)|Rp!def>;
for k in [1..#Pure9Indices] do
 i:=Pure9Indices[k];xyco,n,r2,jco,nr3,orbitLength:=Explode(records[i]);xy:=[Evaluate(Rp!v,a9):v in xyco];x,y:=Explode(xy);
 r:=(x^2*y-x*y+y-1)/(x*(x*y-1));s:=(x*y-y+1)/(x*y);b:=r*s*(r-1);c:=s*(r-1);E:=EllipticCurve([1-c,-b,-b,K9!0,K9!0]);P:=E![0,0];assert TateNormalForm_xy(E,P) eq xy;
 aa,ee:=PairClass(xy,9);assert aa eq h5(L!Pure9Preimages[k]);ss:=SecondClasses(E,P,9);assert #ss eq r2;
 for v in ss do
  pass,im,pre:=ClassCondition(v[1]);fprintf out,"[%o,%o,%o,%o,%o,%o,%o,%o],\n",[1,i,v[3]],[9],[v[2]],Eltseq(v[1]),v[1] in H5 select 1 else 0,pass select 1 else 0,im,pre;Flush(out);
 end for;
 print "NONIC54_SECOND_PURE9",i,#ss,Cputime(tp);
end for;
fprintf out,"]\ncomplete=True\n";delete out;
printf "NONIC54_SECOND_SURVIVOR_TEST_DONE lower=3 pure9=%o complete=true time=%o\n",#Pure9Indices,Cputime(tp);quit;
