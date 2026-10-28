SetColumns(0);SetQuitOnError(true);SetMemoryLimit(32*10^9);tt:=Cputime();
N:=87;gg:=33;prec:=190;h:=16;
assert (h^7 mod N) eq 1 and h ne 1;
MS:=CuspidalSubspace(ModularSymbolsH(N,[h],2,1));assert Dimension(MS) eq gg;
printf "NONIC87_QUOTIENT_SYMBOLS dimension33=true H_index4=true H_order7=true time=%o\n",Cputime(tt);
decomp:=NewformDecomposition(MS);printf "NONIC87_QUOTIENT_CONSTITUENTS dims=%o time=%o\n",[Dimension(A):A in decomp],Cputime(tt);
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
PS<q>:=PowerSeriesRing(Rationals(),prec);Q0:=[];
for A in decomp do
 lev:=Level(AssociatedNewSpace(A));assert lev in [29,87];startRow:=#Q0+1;printf "NONIC87_CONSTITUENT_NEWLEVEL level=%o\n",lev;ff:=qEigenform(A,prec);KF:=BaseRing(Parent(ff));
 if KF cmpeq Rationals() then Append(~Q0,PS!ff);else
  dd:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(ff,n)):n in [0..prec-1]];assert &and[#z eq dd:z in cofs];
  printf "NONIC87_QUOTIENT_EIGENFORM constituent_dim=%o coefficient_dim=%o time=%o\n",Dimension(A),dd,Cputime(tt);
  for k in [1..dd] do Append(~Q0,PS![z[k]:z in cofs]+O(q^prec));end for;
 end if;
if lev eq 29 then oldRows:=Q0[startRow..#Q0];for ff0 in oldRows do Append(~Q0,PS![n mod 3 eq 0 select Coefficient(ff0,n div 3) else 0:n in [0..prec-1]]+O(q^prec));end for;end if;
end for;assert #Q0 eq gg;
QM0:=Matrix(Rationals(),[[Coefficient(f,n):n in [0..prec-1]]:f in Q0]);assert Rank(QM0) eq gg;
den:=LCM([Denominator(a):a in Eltseq(QM0)]);QM:=ChangeRing(Saturation(ChangeRing(den*QM0,Integers())),Rationals());
assert Nrows(QM) eq gg and RowSpace(QM) eq RowSpace(QM0) and Rank(ChangeRing(QM,GF(2))) eq gg;
Q:=[PS!Eltseq(v)+O(q^prec):v in Rows(QM)];pairs:=[];for i in [1..gg],j in [i..gg] do Append(~pairs,<i,j>);end for;
A:=Matrix(Rationals(),[[Coefficient(Q[z[1]]*Q[z[2]],n):n in [0..prec-1]]:z in pairs]);assert Rank(A) eq 3*gg-3;
A2:=ChangeRing(A,GF(2));assert Rank(A2) eq 3*gg-3;rels:=BasisMatrix(Kernel(A2));assert Nrows(rels) eq 465;
out:="C87/data/quotient33/canonical_data.m";
PrintFile(out,Sprintf("N:=87;genus:=33;H_index:=4;precision:=190;qcoeffs:=Matrix(Integers(),33,190,%o);pairs:=%o;QuadRelations:=Matrix(GF(2),465,561,%o);",Eltseq(ChangeRing(QM,Integers())),pairs,Eltseq(rels)):Overwrite:=true);
P:=ProjectiveSpace(GF(2),gg-1);R:=CoordinateRing(P);mons:=[R.(z[1])*R.(z[2]):z in pairs];eqs:=[&+[v[i]*mons[i]:i in [1..#pairs]]:v in Rows(rels)];C:=Curve(P,eqs);assert Dimension(C) eq 1 and Degree(C) eq 64 and ArithmeticGenus(C) eq 33;
printf "NONIC87_QUOTIENT_CANONICAL_READY genus33=true degree64=true relations465=true time=%o\n",Cputime(tt);quit;
