// Prepared only: start through the owned finite-guard wrapper after root allocation.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
N:=71;g:=26;h:=Integers()!((GF(N)!PrimitiveRoot(N))^5);
M:=ModularSymbolsH(N,[h],2,1);MS:=CuspidalSubspace(M);assert Dimension(MS) eq g;
printf "NONIC71_QUOTIENT_SYMBOLS dimension=%o H_index=5 H_order=7 time=%o\n",g,Cputime(tt);
decomp:=NewformDecomposition(MS);printf "NONIC71_QUOTIENT_NEWFORMS dimensions=%o time=%o\n",[Dimension(A):A in decomp],Cputime(tt);
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
PS0<qq>:=PowerSeriesRing(Rationals(),130);Q0:=[];
for A in decomp do
 ff:=qEigenform(A,130);KF:=BaseRing(Parent(ff));
 if KF cmpeq Rationals() then Append(~Q0,PS0!ff);else
  dd:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(ff,n)):n in [0..129]];assert &and[#z eq dd:z in cofs];
  printf "NONIC71_QUOTIENT_EIGENFORM dimension=%o coefficient_dimension=%o time=%o\n",Dimension(A),dd,Cputime(tt);
  for k in [1..dd] do Append(~Q0,PS0![z[k]:z in cofs]+O(qq^130));end for;
 end if;
end for;assert #Q0 eq g;
QM0:=Matrix(Rationals(),[[Coefficient(f,n):n in [0..129]]:f in Q0]);den:=LCM([Denominator(a):a in Eltseq(QM0)]);
QM:=ChangeRing(Saturation(ChangeRing(den*QM0,Integers())),Rationals());assert Nrows(QM) eq g and RowSpace(QM) eq RowSpace(QM0) and Rank(ChangeRing(QM,GF(2))) eq g;
PS<q>:=PowerSeriesRing(Rationals(),130);Q:=[PS!Eltseq(v)+O(q^130):v in Rows(QM)];
pairs:=[];for i in [1..g] do for j in [i..g] do Append(~pairs,<i,j>);end for;end for;assert #pairs eq 351;
A:=Matrix(Rationals(),[[Coefficient(Q[z[1]]*Q[z[2]],n):n in [0..129]]:z in pairs]);
assert Rank(A) eq 75;A2:=ChangeRing(A,GF(2));assert Rank(A2) eq 75;relations:=BasisMatrix(Kernel(A2));assert Nrows(relations) eq 276;
f:="P71/data/quotient26/canonical_data.m";
PrintFile(f,Sprintf("N:=71;H_index:=5;genus:=26;precision:=130;qcoeffs:=Matrix(Integers(),26,130,%o);pairs:=%o;QuadRelations:=Matrix(GF(2),276,351,%o);",Eltseq(ChangeRing(QM,Integers())),pairs,Eltseq(relations)):Overwrite:=true);
printf "NONIC71_QUOTIENT_QUADRICS rank=75 independent_relations=276 time=%o\n",Cputime(tt);
P:=ProjectiveSpace(GF(2),g-1);R:=CoordinateRing(P);mons:=[R.(z[1])*R.(z[2]):z in pairs];eqs:=[&+[v[i]*mons[i]:i in [1..351]]:v in Rows(relations)];
C:=Curve(P,eqs);assert Dimension(C) eq 1 and Degree(C) eq 50 and ArithmeticGenus(C) eq g;
printf "NONIC71_QUOTIENT_CANONICAL_CURVE dimension=1 degree=50 arithmetic_genus=26 time=%o\n",Cputime(tt);
PrintFile(f,"print \"NONIC71_QUOTIENT_CANONICAL_DATA_READY\";");print "NONIC71_QUOTIENT_CANONICAL_PROBE_COMPLETE";quit;
