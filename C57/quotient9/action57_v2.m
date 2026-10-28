SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
N:=57;h:=4;gg:=9;pp:=2;prec:=1100;small:=120;
assert h^9 mod N eq 1 and h ne 1;
MS:=CuspidalSubspace(ModularSymbolsH(N,[h],2,1));assert Dimension(MS) eq gg;
Dcmp:=NewformDecomposition(MS);
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
rows:=[];constits:=[];signs:=[];minusSpace:=Dcmp[1];
for A in Dcmp do
 lev:=Level(AssociatedNewSpace(A));assert N mod lev eq 0;
 ff:=qEigenform(A,prec);KF:=BaseRing(Parent(ff));block:=[];ch:=(Coefficient(ff,5)^2-Coefficient(ff,25))/5;assert ch in [KF|1,-1];sgn:=Rationals()!ch;if sgn eq -1 then assert Dimension(A) eq 4 and LRatio(A,1) ne 0;minusSpace:=A;end if;
 if KF cmpeq Rationals() then Append(~block,[Rationals()!Coefficient(ff,n):n in [0..prec-1]]);else
  dd:=#Flat(KF!1);co:=[Flat(KF!Coefficient(ff,n)):n in [0..prec-1]];assert &and[#v eq dd:v in co];
  for k in [1..dd] do Append(~block,[v[k]:v in co]);end for;
 end if;
 for d in Divisors(N div lev) do
  for r in block do Append(~signs,sgn);Append(~rows,[n mod d eq 0 select r[n div d+1] else Rationals()!0:n in [0..prec-1]]);end for;
 end for;
 Append(~constits,<lev,#block,Divisors(N div lev),Dimension(A)>);
 printf "NONIC57_CONSTITUENT level=%o coefficient_dimension=%o old_divisors=%o symbol_dimension=%o time=%o\n",lev,#block,Divisors(N div lev),Dimension(A),Cputime(tt);
end for;assert #rows eq gg;
Qraw:=Matrix(Rationals(),rows);Q0:=Submatrix(Qraw,1,1,gg,small);assert Rank(Q0) eq gg;
den:=LCM([Denominator(a):a in Eltseq(Q0)]);QS:=ChangeRing(Saturation(ChangeRing(den*Q0,Integers())),Rationals());assert Nrows(QS) eq gg and RowSpace(QS) eq RowSpace(Q0);
Tr:=Solution(Q0,QS);full:=Tr*Qraw;assert Submatrix(full,1,1,gg,small) eq QS;
assert &and[Denominator(a) mod pp ne 0:a in Eltseq(full)];F5:=GF(pp);Q5:=ChangeRing(full,F5);assert Rank(Submatrix(Q5,1,1,gg,small)) eq gg;
E,U:=EchelonForm(Submatrix(Q5,1,1,gg,small));assert U*Submatrix(Q5,1,1,gg,small) eq E;
Q5:=U*Q5;
PS<q>:=PowerSeriesRing(F5,prec);Q:=[PS!Eltseq(v)+O(q^prec):v in Rows(Q5)];
pairs:=[];for i in [1..gg],j in [i..gg] do Append(~pairs,<i,j>);end for;
A:=Matrix(F5,[[Coefficient(Q[z[1]]*Q[z[2]],n):n in [0..small-1]]:z in pairs]);assert Rank(A) eq 3*gg-3;rels:=BasisMatrix(Kernel(A));assert Nrows(rels) eq 21;
out:="C57/data/quotient9/canonical_data.m";
PrintFile(out,Sprintf("N:=57;genus:=9;p:=2;H_generator:=4;precision:=%o;SmallPrecision:=%o;Constituents:=%o;qcoeffs:=Matrix(GF(2),9,%o,%o);pairs:=%o;QuadRelations:=Matrix(GF(2),21,45,%o);",prec,small,constits,prec,Eltseq(Q5),pairs,Eltseq(rels)):Overwrite:=true);
DS:=Tr*DiagonalMatrix(Rationals(),signs)*Tr^-1;assert &and[Denominator(x) mod 2 ne 0:x in Eltseq(DS)];D5:=U*ChangeRing(DS,F5)*U^-1;assert D5^2 eq 1;
PrintFile("C57/data/quotient9/diamond5_data.m",Sprintf("N:=57;genus:=9;p:=2;Diamond5Matrix:=Matrix(GF(2),9,9,%o);",Eltseq(D5)):Overwrite:=true);
printf "NONIC57_DIAMOND5_EXACT canonical_rank=%o rational_trace=%o time=%o\n",Rank(D5-IdentityMatrix(F5,gg)),Trace(DS),Cputime(tt);
assert Dimension(minusSpace) eq 4;ords:=[];I4:=IdentityMatrix(Rationals(),4);
for ell in [2,5,7,11,13,17,23,29,31,37,41,43]do
 TL:=HeckeOperator(minusSpace,ell);DL:=(TL^2-HeckeOperator(minusSpace,ell^2))/ell;vv:=Integers()!Determinant(I4-TL+ell*DL);Append(~ords,<ell,vv>);printf "NONIC57_PRYM_ORDER prime=%o order=%o gcd=%o time=%o\n",ell,vv,GCD([r[2]:r in ords]),Cputime(tt);
end for;
PrintFile("C57/data/quotient9/prym_order_bounds.m",Sprintf("N:=57;PrymDimension:=4;PrymOrders:=%o;PrymRationalOrderBound:=%o;",ords,GCD([r[2]:r in ords])):Overwrite:=true);
print "NONIC57_DIAMOND_AND_PRYM_COMPLETE";quit;
