SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);SetSeed(20261027);tt:=Cputime();F:=GF(2);g:=252;h:=588;
places:=eval Read("C57/data/source85/rr_places_data.txt");cases:=eval Read("C57/data/source85/rr_correct_divisors_data.txt");ss:=Split(Read("C57/data/source85/ring24_tensor.txt"));hdr:=Split(ss[1]," ");assert StringToInteger(hdr[1])eq 31878 and StringToInteger(hdr[2])eq 588;
rows:=[];for k in[2..#ss]do vv:=Intseq(StringToInteger(ss[k],2),2);vv cat:=[0:i in[1..588-#vv]];Append(~rows,Reverse(vv));end for;Tensor:=Matrix(F,rows);assert Nrows(Tensor)eq 31878 and Ncols(Tensor)eq 588;
IX:=Matrix(Integers(),g,g,[0:i in[1..g*g]]);for i,j in[1..g]do a:=Min(i,j);b:=Max(i,j);IX[i,j]:=(a-1)*253-((a-1)*a)div 2+b-a+1;end for;
function Multiply(v)
 M:=ZeroMatrix(F,g,h);for i in Support(v),j in[1..g]do M[j]+:=Tensor[IX[i,j]];end for;return M;
end function;
function Conditions(ts,sign)
 zz:=[z:z in ts|sign*z[2]gt 0];data:=[[F|c:c in &cat[&cat[places[z[1]][3][j][n]:n in[1..Abs(z[2])]]:z in zz]]:j in[1..g]];return Matrix(F,data);
end function;
results:=[* *];
for entry in cases do
 mode,ind,ts:=Explode(entry);plus:=Conditions(ts,1);minus:=Conditions(ts,-1);degree:=&+[Integers()|places[z[1]][1]*z[2]:z in ts|z[2]gt 0];assert Rank(plus)eq degree and Rank(minus)eq degree;V:=BasisMatrix(Kernel(plus));assert Nrows(V)eq g-degree;
 ff:=Random(VectorSpace(F,Nrows(V)))*V;while ff eq 0 do ff:=Random(VectorSpace(F,Nrows(V)))*V;end while;MF:=Multiply(ff);assert Rank(MF)eq g;NK:=Transpose(BasisMatrix(Kernel(Transpose(MF))));assert Ncols(NK)eq h-g and MF*NK eq 0;
 test:=minus;used:=[];hhlist:=[];
 for trial in[1..8]do
  hh:=Random(VectorSpace(F,Nrows(V)))*V;assert hh*plus eq 0;MH:=Multiply(hh);test:=HorizontalJoin(test,MH*NK);Append(~hhlist,Eltseq(hh));Append(~used,Rank(test));if Rank(test)eq g then break;end if;
 end for;
 assert ff*plus eq 0;Append(~results,[*mode,ind,ts,degree,Rank(plus),Rank(minus),Eltseq(ff),hhlist,used,Rank(MF),Rank(test)*]);
 printf "NONIC57_SOURCE85_RR24_CORRECT_CASE case=%o index=%o degree=%o positive_rank=%o negative_rank=%o multiplication_rank252=true final_rank=%o H_trials=%o time=%o\n",mode,ind,degree,Rank(plus),Rank(minus),Rank(test),used,Cputime(tt);
end for;
Write("C57/data/source85/rr24_correct_results.txt",Sprint(results):Overwrite:=true);assert &and[v[11]eq g:v in results];printf "NONIC57_SOURCE85_RR24_CORRECT_ALL_EXCLUSIONS_COMPLETE cases=%o time=%o\n",#results,Cputime(tt);quit;
