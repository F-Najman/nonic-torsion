load "C77/data/retained_degree_independent_kernels.m";
rt:=Cputime();n:=37;D:=9;V:=VectorSpace(GF(2),n*n);eqs:=[];freshBM:=BasisMatrix(U77);
i0:=Index(labels77,CuspLabels[SourceCusp]);assert i0 gt 0;
WW,wwmap:=DifferentialSpace(8*Divisor(cs77[i0]));Wnew8:=sub<VectorSpace(GF(2),37)|[Eltseq(wwmap(v)@@om77):v in Basis(WW)]>;
Anew8:=Transpose(BasisMatrix(Nullspace(Transpose(BasisMatrix(Wnew8)))));assert Ncols(Anew8) eq 8;
for pair in [<CartierMatrix,CM77>,<MatrixC2,M2>,<MatrixC3,M3>] do
 Cold:=pair[1];Cnew:=pair[2];
 for i,j in [1..n] do
  row:=V!0;for k in [1..n] do row[(k-1)*n+j]+:=Cold[i,k];row[(i-1)*n+k]+:=Cnew[k,j];end for;Append(~eqs,row);
 end for;
end for;
Wold8:=sourcekernels[8];
for r in [1..Nrows(Wold8)],c in [1..8] do
 row:=V!0;for i,j in [1..n] do row[(i-1)*n+j]:=Wold8[r,i]*Anew8[j,c];end for;Append(~eqs,row);
end for;
Ktrans:=Nullspace(Transpose(Matrix(eqs)));assert Dimension(Ktrans) eq 1;
Ttrans:=Matrix(GF(2),n,n,Eltseq(Basis(Ktrans)[1]));assert Rank(Ttrans) eq 37;
assert CartierMatrix*Ttrans eq Ttrans*CM77 and MatrixC2*Ttrans eq Ttrans*M2 and MatrixC3*Ttrans eq Ttrans*M3;
assert RowSpace(Wold8*Ttrans) eq Wnew8 and RowSpace(ImageBasis*Ttrans) eq U77;
printf "NONIC77_BASIS_ALIGNMENT_PASS intertwiner_dimension=1 invertible=true anchor_m=8 time=%o\n",Cputime(rt);
csRat:=[cs77[Index(labels77,l)]:l in CuspLabels];
csCub:=[cys[i]:i in [1..#cys]|cywidths[i] eq 11 and Degree(cys[i]) eq 3];assert #csCub eq 2;
fq,gq:=TwoGenerators(csCub[1]);assert Valuation(deck2(fq),csCub[2]) gt 0 and Valuation(deck2(gq),csCub[2]) gt 0;
VO:=VectorSpace(GF(2),37);blocks:=AssociativeArray();Ablocks:=AssociativeArray();Wblocks:=AssociativeArray();
for i in [1..6] do
 for m in [1..8] do
  AA:=Ttrans^-1*ColumnSubmatrix(basisjets[i],m);assert Nullspace(AA) eq RowSpace(kernelspaces[<i,m>]*Ttrans);
  Ablocks[<i,m>]:=AA;Wblocks[<i,m>]:=BasisMatrix(Nullspace(AA));blocks[<i,m>]:=freshBM*AA;
 end for;
end for;
WW,wwmap:=DifferentialSpace(9*Divisor(cs77[i0]));W9:=sub<VO|[Eltseq(wwmap(v)@@om77):v in Basis(WW)]>;
AA9:=Transpose(BasisMatrix(Nullspace(Transpose(BasisMatrix(W9)))));while Ncols(AA9) lt 9 do AA9:=HorizontalJoin(AA9,ZeroMatrix(GF(2),37,1));end while;assert Ncols(AA9) eq 9 and Nullspace(AA9) eq W9;
for i in [1..6] do
 matches:=[];
 for e2 in [0,1],e3 in [0,1,2] do
  k:=SourceCusp;for r in [1..e2] do k:=CuspPermutation2[k];end for;for r in [1..e3] do k:=CuspPermutation3[k];end for;
  if k eq i then Append(~matches,<e2,e3>);end if;
 end for;
 assert #matches eq 1;ee:=matches[1];MM:=M2^ee[1]*M3^ee[2];AA:=MM*AA9;
 Ablocks[<i,9>]:=AA;Wblocks[<i,9>]:=BasisMatrix(Nullspace(AA));blocks[<i,9>]:=freshBM*AA;
end for;
printf "NONIC77_NEW_RATIONAL_M9_PASS rank=%o time=%o\n",Rank(freshBM*AA9),Cputime(rt);
for m in [1..3] do
 WW,wwmap:=DifferentialSpace(m*Divisor(csCub[1]));W:=sub<VO|[Eltseq(wwmap(v)@@om77):v in Basis(WW)]>;
 AA:=Transpose(BasisMatrix(Nullspace(Transpose(BasisMatrix(W)))));while Ncols(AA) lt 3*m do AA:=HorizontalJoin(AA,ZeroMatrix(GF(2),37,1));end while;assert Ncols(AA) eq 3*m and Nullspace(AA) eq W;
 for k in [1,2] do
  MAT:=M2^(k-1)*AA;idx:=6+k;Ablocks[<idx,m>]:=MAT;Wblocks[<idx,m>]:=BasisMatrix(Nullspace(MAT));blocks[<idx,m>]:=freshBM*MAT;
 end for;
 printf "NONIC77_NEW_CUBIC_RR m=%o dimension=%o restriction_rank=%o time=%o\n",m,Dimension(W),Rank(freshBM*AA),Cputime(rt);
end for;
function Weighted77(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else [];end if;
 ans:=[];for m in [0..d div ds[1]] do for vv in Weighted77(ds[2..#ds],d-m*ds[1]) do Append(~ans,[m] cat vv);end for;end for;return ans;
end function;
ds9:=[1,1,1,1,1,1,3,3];vectors:=Weighted77(ds9,9);assert #vectors eq 3098;bad:=[];
for idx in [1..#vectors] do
 vv:=vectors[idx];AA:=ZeroMatrix(GF(2),32,0);
 for i in [1..8] do if vv[i] gt 0 then AA:=HorizontalJoin(AA,blocks[<i,vv[i]>]);end if;end for;
 assert Ncols(AA) eq 9;if Rank(AA) ne 9 then Append(~bad,vv);print "NONIC77_FI9_EXCEPTION",idx,Rank(AA),vv;end if;
end for;
if not assigned OutFile then OutFile:="C77/data/fi9_data.json";end if;
out:=Open(OutFile,"w");fprintf out,"{\"genus\":37,\"degree\":9,\"Cartier\":%o,\"U_dimension\":32,\"U_basis\":%o,\"T_old_to_new\":%o,\"deck2\":%o,\"deck3\":%o,\"anchor_kernel\":%o,\"anchor_restriction\":%o,\"degrees\":%o,\"blocks\":[",[Integers()!x:x in Eltseq(CM77)],[Integers()!x:x in Eltseq(freshBM)],[Integers()!x:x in Eltseq(Ttrans)],[Integers()!x:x in Eltseq(M2)],[Integers()!x:x in Eltseq(M3)],[Integers()!x:x in Eltseq(BasisMatrix(Wnew8))],[Integers()!x:x in Eltseq(Anew8)],ds9;
first:=true;
for i in [1..8] do
 for m in [1..(9 div ds9[i])] do
  if not first then fprintf out,",";end if;first:=false;AA:=Ablocks[<i,m>];WW:=Wblocks[<i,m>];
  fprintf out,"{\"cusp\":%o,\"multiplicity\":%o,\"kernel_dimension\":%o,\"kernel\":%o,\"restriction\":%o}",i,m,Nrows(WW),[Integers()!x:x in Eltseq(WW)],[Integers()!x:x in Eltseq(AA)];
 end for;
end for;
fprintf out,"],\"tested\":3098,\"exceptions\":%o}\n",bad;delete out;
printf "NONIC77_FI9_COMPLETE tested=3098 exceptions=%o complete=true time=%o\n",#bad,Cputime(rt);
quit;
