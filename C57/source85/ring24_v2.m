SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
load "C57/data/source85/canonical_unit_data.m";
g:=85;prec:=1100;F:=GF(2);PS<q>:=PowerSeriesRing(F,prec);Q:=[PS!Eltseq(v)+O(q^prec):v in Rows(qcoeffs)];p1:=[];for i in[1..g]do for j in[i..g]do Append(~p1,<i,j>);end for;end for;prod:=[Q[z[1]]*Q[z[2]]:z in p1];A:=Matrix(F,[[Coefficient(v,n):n in[0..480]]:v in prod]);assert Rank(A)eq 252;E:=EchelonForm(Transpose(A));idx:=[Index(Eltseq(E[i]),F!1):i in[1..252]];assert #Seqset(idx)eq 252;B:=[prod[i]:i in idx];QB:=Matrix(F,[[Coefficient(v,n):n in[0..prec-1]]:v in B]);C12:=Solution(Submatrix(QB,1,1,252,481),A);All12:=Matrix(F,[[Coefficient(v,n):n in[0..prec-1]]:v in prod]);assert C12*QB eq All12;
printf "NONIC57_SOURCE85_RING2_COMPLETE dimension252=true all3655products=true time=%o\n",Cputime(tt);
p2:=[];for i in[1..252]do for j in[i..252]do Append(~p2,<i,j>);end for;end for;prod4:=[B[z[1]]*B[z[2]]:z in p2];A4:=Matrix(F,[[Coefficient(v,n):n in[0..960]]:v in prod4]);assert Rank(A4)eq 588;E4:=EchelonForm(Transpose(A4));idx4:=[Index(Eltseq(E4[i]),F!1):i in[1..588]];assert #Seqset(idx4)eq 588;B4:=[prod4[i]:i in idx4];QB4:=Matrix(F,[[Coefficient(v,n):n in[0..prec-1]]:v in B4]);C24:=Solution(Submatrix(QB4,1,1,588,961),A4);All24:=Matrix(F,[[Coefficient(v,n):n in[0..prec-1]]:v in prod4]);assert C24*QB4 eq All24;
printf "NONIC57_SOURCE85_RING4_COMPLETE dimension588=true all31878products=true time=%o\n",Cputime(tt);
LS<z>:=LaurentSeriesRing(F,prec);xx:=z^-1+&+[LS|F!UnitXSeries[n+1]*z^n:n in[0..#UnitXSeries-1]]+O(z^prec);yy:=z^-1+&+[LS|F!UnitYSeries[n+1]*z^n:n in[0..#UnitYSeries-1]]+O(z^prec);LL:=[LS!Eltseq(v)+O(z^prec):v in Rows(QB)];
function URel(uu)
 M:=Matrix(F,[[Coefficient(v,n):n in[0..1098]]:v in LL cat[uu*v:v in LL]]);return BasisMatrix(Kernel(M));
end function;
RX:=URel(xx);RY:=URel(yy);assert Nrows(RX)eq 215 and Nrows(RY)eq 213;
PrintFile("C57/data/source85/ring24_data.m",Sprintf("N:=57;genus:=85;precision:=1100;Ring1Pairs:=%o;Ring2MonomialIndices:=%o;Ring2Pairs:=%o;Ring4MonomialIndices:=%o;Ring2qCoefficients:=Matrix(GF(2),252,1100,%o);Ring4qCoefficients:=Matrix(GF(2),588,1100,%o);Ring2UnitXRelations:=Matrix(GF(2),215,504,%o);Ring2UnitYRelations:=Matrix(GF(2),213,504,%o);",p1,idx,p2,idx4,Eltseq(QB),Eltseq(QB4),Eltseq(RX),Eltseq(RY)):Overwrite:=true);
fn:="C57/data/source85/ring12_tensor.txt";PrintFile(fn,"3655 252":Overwrite:=true);for v in Rows(C12)do PrintFile(fn,&cat[IntegerToString(Integers()!c):c in Eltseq(v)]);end for;
fn:="C57/data/source85/ring24_tensor.txt";PrintFile(fn,"31878 588":Overwrite:=true);for v in Rows(C24)do PrintFile(fn,&cat[IntegerToString(Integers()!c):c in Eltseq(v)]);end for;
printf "NONIC57_SOURCE85_RING24_TENSOR_COMPLETE R2unitrelations=215,213=true time=%o\n",Cputime(tt);quit;
