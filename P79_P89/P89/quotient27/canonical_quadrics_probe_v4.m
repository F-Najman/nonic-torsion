SetColumns(0);SetQuitOnError(true);SetMemoryLimit(32*10^9);tt:=Cputime();
N:=89;h:=Integers()!((GF(89)!PrimitiveRoot(89))^4);
M:=ModularSymbolsH(89,[h],2,1);MS:=CuspidalSubspace(M);assert Dimension(MS) eq 27;
printf "NONIC89_QUOTIENT_SYMBOLS dimension=27 H_index=4 H_order=11 time=%o\n",Cputime(tt);
decomp:=NewformDecomposition(MS);printf "NONIC89_QUOTIENT_NEWFORMS dimensions=%o time=%o\n",[Dimension(A):A in decomp],Cputime(tt);
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
PS0<qq>:=PowerSeriesRing(Rationals(),130);Q0:=[];
for A in decomp do
 ff:=qEigenform(A,130);KF:=BaseRing(Parent(ff));
 if KF cmpeq Rationals() then Append(~Q0,PS0!ff);else
  dd:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(ff,n)):n in [0..129]];assert &and[#z eq dd:z in cofs];printf "NONIC89_QUOTIENT_EIGENFORM dimension=%o absolute_coefficient_dimension=%o coefficient_ring_type=%o time=%o\n",Dimension(A),dd,Type(KF),Cputime(tt);
  for k in [1..dd] do Append(~Q0,PS0![z[k]:z in cofs]+O(qq^130));end for;
 end if;
end for;assert #Q0 eq 27;
QM0:=Matrix(Rationals(),[[Coefficient(f,n):n in [0..129]]:f in Q0]);den:=LCM([Denominator(a):a in Eltseq(QM0)]);
QM:=ChangeRing(Saturation(ChangeRing(den*QM0,Integers())),Rationals());assert Nrows(QM) eq 27 and RowSpace(QM) eq RowSpace(QM0) and Rank(ChangeRing(QM,GF(2))) eq 27;
PS<q>:=PowerSeriesRing(Rationals(),130);Q:=[PS!Eltseq(v)+O(q^130):v in Rows(QM)];
PrintFile("P79_P89/data/P89/quotient27/q_basis_data.m",Sprintf("N:=89;H_index:=4;genus:=27;precision:=130;qcoeffs:=Matrix(Integers(),27,130,%o);",Eltseq(ChangeRing(QM,Integers()))):Overwrite:=true);
printf "NONIC89_QUOTIENT_Q_BASIS integral=true mod2_rank=27 precision=130 time=%o\n",Cputime(tt);
pairs:=[];for i in [1..27] do for j in [i..27] do Append(~pairs,<i,j>);end for;end for;assert #pairs eq 378;
A:=Matrix(Rationals(),[[Coefficient(Q[z[1]]*Q[z[2]],n):n in [0..129]]:z in pairs]);
assert Rank(A) eq 78;A2:=ChangeRing(A,GF(2));assert Rank(A2) eq 78;relations:=BasisMatrix(Kernel(A2));assert Nrows(relations) eq 300;
f:="P79_P89/data/P89/quotient27/canonical_data.m";
PrintFile(f,Sprintf("N:=89;H_index:=4;genus:=27;precision:=130;qcoeffs:=Matrix(Integers(),27,130,%o);pairs:=%o;QuadRelations:=Matrix(GF(2),300,378,%o);",Eltseq(ChangeRing(QM,Integers())),pairs,Eltseq(relations)):Overwrite:=true);
printf "NONIC89_QUOTIENT_QUADRICS rank=78 independent_relations=300 time=%o\n",Cputime(tt);
P:=ProjectiveSpace(GF(2),26);R:=CoordinateRing(P);mons:=[R.(z[1])*R.(z[2]):z in pairs];eqs:=[&+[v[i]*mons[i]:i in [1..378]]:v in Rows(relations)];
C:=Curve(P,eqs);assert Dimension(C) eq 1;
printf "NONIC89_QUOTIENT_CANONICAL_CURVE dimension=%o degree=%o arithmetic_genus=%o time=%o\n",Dimension(C),Degree(C),ArithmeticGenus(C),Cputime(tt);
assert Degree(C) eq 52 and ArithmeticGenus(C) eq 27;
PrintFile(f,"print \"NONIC89_QUOTIENT_CANONICAL_DATA_READY\";");
print "NONIC89_QUOTIENT_CANONICAL_QUADRICS_PROBE_COMPLETE";quit;
