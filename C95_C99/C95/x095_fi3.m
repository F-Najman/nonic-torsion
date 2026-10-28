SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
CQ:=SmallModularCurve(19);jraw:=jInvariant(CQ,19);CC:=ChangeRing(CQ,GF(2));assert IsNonsingular(CC) and Genus(CC) eq 1;
RC:=CoordinateRing(Ambient(CC));jn:=RC!Numerator(jraw);jd:=RC!Denominator(jraw);
K,kmap:=AlgorithmicFunctionField(FunctionField(CC));jK:=kmap(FunctionField(CC)!(jn/jd));
Fv<v>:=FunctionField(GF(2));j5:=v^6/(v+1);j5dual:=v^6/(v+1)^5;assert Evaluate(ClassicalModularPolynomial(5),[j5,j5dual]) eq 0;
R<x>:=PolynomialRing(K);f:=x^6+jK*x+jK;assert IsIrreducible(f);AF<u>:=FunctionField(f);assert Genus(AF) eq 9;
cs:=Poles(AF!jK);widths:=[-Valuation(AF!jK,c):c in cs];assert #cs eq 4 and &and[Degree(c) eq 1:c in cs] and Seqset(widths) eq {1,5,19,95};
printf "NONIC_95_X0_MODEL genus9=true cusp_widths=%o time=%o\n",widths,Cputime(tt);
Om,om:=SpaceOfHolomorphicDifferentials(AF);bs:=[om(v):v in Basis(Om)];assert Dimension(Om) eq 9;
Cart:=Matrix(GF(2),[Eltseq(Cartier(w)@@om):w in bs]);U:=RowSpace(Cart-IdentityMatrix(GF(2),9));BM:=BasisMatrix(U);udim:=Dimension(U);
printf "NONIC_95_X0_CARTIER rank=%o image_dimension=%o time=%o\n",Rank(Cart),udim,Cputime(tt);
if not assigned OutFile then OutFile:="C95_C99/data/C95/fi3_data.json";end if;
out:=Open(OutFile,"w");
fprintf out,"{\"genus\":9,\"degree\":3,\"Cartier\":%o,\"U_dimension\":%o,\"U_basis\":%o,\"cusp_widths\":%o,\"kernels\":[",[Integers()!x:x in Eltseq(Cart)],udim,[Integers()!x:x in Eltseq(BM)],widths;
V:=VectorSpace(GF(2),9);blocks:=AssociativeArray();first:=true;
for i in [1..4] do
 for m in [1..3] do
  W,wm:=DifferentialSpace(m*Divisor(cs[i]));WW:=sub<V|[V!Eltseq(wm(v)@@om):v in Basis(W)]>;
  ann:=Nullspace(Transpose(BasisMatrix(WW)));AA:=Transpose(BasisMatrix(ann));assert Ncols(AA) eq m and Nullspace(AA) eq WW;
  blocks[<i,m>]:=BM*AA;
  if not first then fprintf out,",";end if;first:=false;
  fprintf out,"{\"cusp\":%o,\"multiplicity\":%o,\"kernel_dimension\":%o,\"kernel\":%o,\"restriction\":%o}",i,m,Dimension(WW),[Integers()!x:x in Eltseq(BasisMatrix(WW))],[Integers()!x:x in Eltseq(AA)];
 end for;
end for;
fprintf out,"],\"configurations\":[";tested:=0;bad:=[];first:=true;
for m1,m2,m3 in [0..3] do
 m4:=3-m1-m2-m3;if m4 lt 0 then continue;end if;vv:=[m1,m2,m3,m4];AA:=ZeroMatrix(GF(2),udim,0);
 for i in [1..4] do if vv[i] gt 0 then AA:=HorizontalJoin(AA,blocks[<i,vv[i]>]);end if;end for;
 assert Ncols(AA) eq 3;rk:=Rank(AA);tested+:=1;
 if rk ne 3 then Append(~bad,vv);end if;
 if not first then fprintf out,",";end if;first:=false;
 fprintf out,"{\"multiplicities\":%o,\"rank\":%o}",vv,rk;
end for;
fprintf out,"]}\n";delete out;assert tested eq 20;
printf "NONIC_95_X0_FI3_COMPLETE tested=%o exceptions=%o bad=%o time=%o\n",tested,#bad,bad,Cputime(tt);
quit;
