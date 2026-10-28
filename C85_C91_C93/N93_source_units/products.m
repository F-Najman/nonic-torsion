SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
load "C85_C91_C93/data/N93_source_units/canonical_unit_data.m";
P<z>:=PolynomialRing(GF(2));Q:=[P!Eltseq(v):v in Rows(qcoeffs)];g:=genus;pr:=precision;
idx:=N^2*&*[Rationals()|1-1/p^2:p in PrimeDivisors(N)]/2;assert IsIntegral(idx);sturm:=Floor(idx/3);assert pr-1 gt sturm;
pairs:=[<i,j>:i,j in[1..g]|i le j];rows:=[];
for ij in pairs do
 f:=(Q[ij[1]]*Q[ij[2]])mod z^pr;c:=Eltseq(f);c cat:=[GF(2)!0:i in[#c+1..pr]];Append(~rows,c);
end for;
M:=Matrix(GF(2),rows);delete rows;rr:=BasisMatrix(RowSpace(M));assert Nrows(rr)eq 3*g-3;
PrintFile("C85_C91_C93/data/N93_source_units/R2_qbasis_data.m",Sprintf("N:=%o;genus:=%o;precision:=%o;SturmWeight4:=%o;R2qBasis:=Matrix(GF(2),%o,%o,%o);Complete:=true;",N,g,pr,sturm,Nrows(rr),pr,Eltseq(rr)):Overwrite:=true);
printf "SOURCE_PRODUCT_COMPLETE level=%o genus=%o dimR2=%o symmetric_products=%o precision=%o SturmWeight4=%o time=%o\n",N,g,Nrows(rr),#pairs,pr,sturm,Cputime(tt);quit;
