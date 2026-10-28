assert udim ge 3;
cs:=Poles(jaf);
assert {*Degree(q):q in cs*} eq {*1^^6,3^^2*};assert &+[-Valuation(jaf,q)*Degree(q):q in cs] eq 480;
ds:=[Degree(q):q in cs];printf "NONIC133_Z_CUSPS degrees=%o jdegree=480 time=%o\n",ds,Cputime(tt);
out:=Open("C133/data/fi3_data_reproduced.m","w");
fprintf out,"N:=133;D:=3;gx:=35;o:=1;udim:=%o;JetFactor:=1;ds:=%o;As:=[];perms:=[];localjets:=AssociativeArray();basisjets:=AssociativeArray();kernelspaces:=AssociativeArray();\n",udim,ds;
fprintf out,"CartierMatrix:=Matrix(GF(2),35,35,%o);ImageBasis:=Matrix(GF(2),%o,35,%o);\n",Eltseq(CM),udim,Eltseq(BasisMatrix(U));
fprintf out,"DiamondMatrix:=Matrix(GF(2),35,35,%o);\n",Eltseq(DM);
VO:=VectorSpace(GF(2),g);jets:=AssociativeArray();
for i in [1..#cs] do
 cols:=[];annprev:=sub<VO|>;
 for m in [1..3 div ds[i]] do
  VV,vvmap:=DifferentialSpace(m*Divisor(cs[i]));W:=sub<VO|[VO!Eltseq(vvmap(w) @@ om):w in Basis(VV)]>;
  ann:=Nullspace(Transpose(BasisMatrix(W)));assert annprev subset ann;
  for v in Basis(ann) do if v notin annprev then Append(~cols,VO!v);annprev:=sub<VO|Basis(annprev) cat [v]>;end if;end for;
  while #cols lt m*ds[i] do Append(~cols,VO!0);end while;
  A:=Transpose(Matrix(GF(2),[Eltseq(v):v in cols]));assert Nullspace(A) eq W;
  J:=BasisMatrix(U)*A;jets[<i,m>]:=J;
  fprintf out,"kernelspaces[<%o,%o>]:=Matrix(GF(2),%o,35,%o);localjets[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);\n",i,m,Dimension(W),Eltseq(BasisMatrix(W)),i,m,Nrows(J),Ncols(J),Eltseq(J);
  Flush(out);printf "NONIC133_Z_RR cusp=%o degree=%o multiplicity=%o rank=%o time=%o\n",i,ds[i],m,Rank(J),Cputime(tt);
 end for;
 fprintf out,"basisjets[%o]:=Matrix(GF(2),35,%o,%o);\n",i,Ncols(A),Eltseq(A);
end for;
fprintf out,"print \"FI_JET_DATA_READY\";\n";delete out;
function WV(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else [];end if;
 result:=[];
 for m in [0..d div ds[1]] do for v in WV(ds[2..#ds],d-m*ds[1]) do Append(~result,[m] cat v);end for;end for;
 return result;
end function;
vs:=WV(ds,3);assert #vs eq 58;bad:=[];
for v in vs do
 A:=ZeroMatrix(GF(2),udim,0);for i in [1..#cs] do if v[i] gt 0 then A:=HorizontalJoin(A,jets[<i,v[i]>]);end if;end for;
 assert Ncols(A) eq 3;if Rank(A) lt 3 then Append(~bad,v);print "NONIC133_Z_FI_EXCEPTION",v,Rank(A);end if;
end for;
assert #bad eq 0;
printf "NONIC133_Z_FI_DONE degree=3 tested=58 exceptions=%o complete=true time=%o\n",#bad,Cputime(tt);quit;
