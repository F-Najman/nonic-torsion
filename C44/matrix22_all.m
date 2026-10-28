SieveOutput:="C44/data/matrix22_all.raw.jsonseq";
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);SetSeed(20261044);tt:=Cputime();
load "C44/data/multiplication_22_data.m";
load "C44/data/sieve_all_inputs.m";
assert #Tensor22Blocks eq 105 and #Std2Exponents eq 105 and #Std4Exponents eq 245;
PF:=PolynomialRing(GF(3));
function EvalJets22(J)
 dd:=J[1];mm:=J[2];LF<a>:=ext<GF(3)|PF!J[3]>;assert Degree(LF) eq dd;PS<T>:=PowerSeriesRing(LF,mm);pp:=[LF!1];for k in [1..dd-1] do Append(~pp,pp[#pp]*a);end for;zs:=[];
 for row in J[4] do cc:=[];for v in row do assert #v eq dd;Append(~cc,&+[LF!v[k]*pp[k]:k in [1..dd]]);end for;assert #cc eq mm;Append(~zs,PS!cc+O(T^mm));end for;assert #zs eq 36;
 powers:=AssociativeArray();for i in [1..36] do powers[<i,0>]:=PS!1;for k in [1..2] do powers[<i,k>]:=powers[<i,k-1>]*zs[i];end for;end for;rows:=[];
 for ee in Std2Exponents do vv:=&*[powers[<i,ee[i]>]:i in [1..36]|ee[i] gt 0];Append(~rows,&cat[Eltseq(Coefficient(vv,k),GF(3)):k in [0..mm-1]]);end for;
 E2:=Matrix(GF(3),rows);assert Nrows(E2) eq 105 and Ncols(E2) eq dd*mm;return E2;
end function;
out:=Open(SieveOutput,"w");fprintf out,"{\"type\":\"header\",\"level\":44,\"p\":3,\"cases\":%o,\"R2\":105,\"R4\":245,\"criterion\":\"R2_times_R2\"}\n",#JetCases;Flush(out);survivors:=[];
for jj in JetCases do
 k:=jj[1];ET:=ZeroMatrix(GF(3),105,0);EU:=ZeroMatrix(GF(3),105,0);for J in jj[2] do ET:=HorizontalJoin(ET,EvalJets22(J));end for;for J in jj[3] do EU:=HorizontalJoin(EU,EvalJets22(J));end for;assert Ncols(ET) eq 54 and Ncols(EU) eq 54;
 assert Rank(ET) eq 54 and Rank(EU) eq 54;VF:=Nullspace(ET);assert Dimension(VF) eq 51;rank0:=0;attempts:=0;
 while attempts lt 8 and rank0 lt 105 do
  attempts+:=1;F:=Random(VF);while IsZero(F) do F:=Random(VF);end while;H:=Random(VF);while IsZero(H) do H:=Random(VF);end while;assert IsZero(F*ET) and IsZero(H*ET);
  MF:=ZeroMatrix(GF(3),105,245);for i in [1..105] do if F[i] ne 0 then MF+:=F[i]*Tensor22Blocks[i];end if;end for;assert Rank(MF) eq 105;
  QQ:=Transpose(BasisMatrix(Nullspace(Transpose(MF))));assert Nrows(QQ) eq 245 and Ncols(QQ) eq 140 and IsZero(MF*QQ);
  MH:=Matrix(GF(3),[Eltseq(H*Tensor22Blocks[i]):i in [1..105]]);AA:=HorizontalJoin(EU,MH*QQ);assert Nrows(AA) eq 105 and Ncols(AA) eq 194;rank0:=Rank(AA);
 end while;
 if rank0 lt 105 then Append(~survivors,k);end if;
 fprintf out,"{\"type\":\"case\",\"case\":%o,\"attempts\":%o,\"rank\":%o,\"F\":%o,\"H\":%o,\"ET2\":%o,\"EU2\":%o,\"quotient_projection\":%o,\"matrix\":%o}\n",k,attempts,rank0,[Integers()!c:c in Eltseq(F)],[Integers()!c:c in Eltseq(H)],[Integers()!c:c in Eltseq(ET)],[Integers()!c:c in Eltseq(EU)],[Integers()!c:c in Eltseq(QQ)],[Integers()!c:c in Eltseq(AA)];Flush(out);
 printf "NONIC44_CANONICAL22_CASE index=%o rank=%o attempts=%o time=%o\n",k,rank0,attempts,Cputime(tt);
end for;
fprintf out,"{\"type\":\"complete\",\"cases\":%o,\"survivors\":%o}\n",#JetCases,survivors;delete out;
printf "NONIC44_CANONICAL22_COMPLETE cases=%o survivors=%o all_rank105=%o time=%o\n",#JetCases,survivors,#survivors eq 0,Cputime(tt);quit;
