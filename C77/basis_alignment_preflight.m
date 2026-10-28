SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
load "C77/data/retained_degree_independent_kernels.m";
n:=37;V:=VectorSpace(GF(2),n*n);eqs:=[];
for C in [CartierMatrix,MatrixC2,MatrixC3] do
 for i,j in [1..n] do
  row:=V!0;
  for k in [1..n] do row[(k-1)*n+j]+:=C[i,k];row[(i-1)*n+k]+:=C[k,j];end for;
  Append(~eqs,row);
 end for;
end for;
for m in [8,4,2,1] do
 W:=sourcekernels[m];A:=ColumnSubmatrix(basisjets[1],m);
 for r in [1..Nrows(W)],c in [1..m] do
  row:=V!0;
  for i,j in [1..n] do row[(i-1)*n+j]:=W[r,i]*A[j,c];end for;
  Append(~eqs,row);
 end for;
 K:=Nullspace(Transpose(Matrix(eqs)));
 printf "NONIC77_BASIS_STABILIZER flag_m=%o dimension=%o time=%o\n",m,Dimension(K),Cputime(tt);
 if Dimension(K) eq 1 then assert Matrix(GF(2),n,n,Eltseq(Basis(K)[1])) eq IdentityMatrix(GF(2),n);break;end if;
end for;
print "NONIC77_BASIS_PREFLIGHT_COMPLETE";quit;
