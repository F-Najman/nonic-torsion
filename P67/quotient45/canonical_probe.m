// Root allocated one small C67 pipeline slot, 8GB, nice10, finite CPU/wall guards.
// Prepared only: start through the owned finite-guard wrapper after root allocation.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
N:=67;g:=45;h:=Integers()!((GF(N)!PrimitiveRoot(N))^11);
M:=ModularSymbolsH(N,[h],2,1);MS:=CuspidalSubspace(M);assert Dimension(MS) eq g;
printf "NONIC67_QUOTIENT_SYMBOLS dimension=%o H_index=11 H_order=3 time=%o\n",g,Cputime(tt);
decomp:=NewformDecomposition(MS);printf "NONIC67_QUOTIENT_NEWFORMS dimensions=%o time=%o\n",[Dimension(A):A in decomp],Cputime(tt);
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
PS0<qq>:=PowerSeriesRing(Rationals(),320);Q0:=[];
for A in decomp do
 ff:=qEigenform(A,320);KF:=BaseRing(Parent(ff));
 if KF cmpeq Rationals() then Append(~Q0,PS0!ff);else
  dd:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(ff,n)):n in [0..319]];assert &and[#z eq dd:z in cofs];
  printf "NONIC67_QUOTIENT_EIGENFORM dimension=%o coefficient_dimension=%o time=%o\n",Dimension(A),dd,Cputime(tt);
  for k in [1..dd] do Append(~Q0,PS0![z[k]:z in cofs]+O(qq^320));end for;
 end if;
end for;assert #Q0 eq g;
QM0:=Matrix(Rationals(),[[Coefficient(f,n):n in [0..319]]:f in Q0]);den:=LCM([Denominator(a):a in Eltseq(QM0)]);
QM:=ChangeRing(Saturation(ChangeRing(den*QM0,Integers())),Rationals());assert Nrows(QM) eq g and RowSpace(QM) eq RowSpace(QM0) and Rank(ChangeRing(QM,GF(2))) eq g;
PS<q>:=PowerSeriesRing(Rationals(),320);Q:=[PS!Eltseq(v)+O(q^320):v in Rows(QM)];
pairs:=[];for i in [1..g] do for j in [i..g] do Append(~pairs,<i,j>);end for;end for;assert #pairs eq 1035;
A:=Matrix(Rationals(),[[Coefficient(Q[z[1]]*Q[z[2]],n):n in [0..319]]:z in pairs]);
assert Rank(A) eq 132;A2:=ChangeRing(A,GF(2));assert Rank(A2) eq 132;relations:=BasisMatrix(Kernel(A2));assert Nrows(relations) eq 903;
f:="P67/data/quotient45/canonical_data.m";
PrintFile(f,Sprintf("N:=67;H_index:=11;genus:=45;precision:=320;qcoeffs:=Matrix(Integers(),45,320,%o);pairs:=%o;QuadRelations:=Matrix(GF(2),903,1035,%o);",Eltseq(ChangeRing(QM,Integers())),pairs,Eltseq(relations)):Overwrite:=true);
printf "NONIC67_QUOTIENT_QUADRICS rank=132 independent_relations=903 time=%o\n",Cputime(tt);
P:=ProjectiveSpace(GF(2),g-1);R:=CoordinateRing(P);mons:=[R.(z[1])*R.(z[2]):z in pairs];eqs:=[&+[v[i]*mons[i]:i in [1..1035]]:v in Rows(relations)];
C:=Curve(P,eqs);assert Dimension(C) eq 1 and Degree(C) eq 88 and ArithmeticGenus(C) eq g;
printf "NONIC67_QUOTIENT_CANONICAL_CURVE dimension=1 degree=88 arithmetic_genus=45 time=%o\n",Cputime(tt);
PrintFile(f,"print \"NONIC67_QUOTIENT_CANONICAL_DATA_READY\";");print "NONIC67_QUOTIENT_CANONICAL_PROBE_COMPLETE";quit;
