SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
for qq in [3,5,7,11,17,23] do
NN:=1;dim:=0;
for chi0 in GaloisConjugacyRepresentatives(FullDirichletGroup(247)) do
 if IsOdd(chi0) or Evaluate(chi0,40) ne 1 or Order(chi0) eq 1 then continue;end if;
 chi:=MinimalBaseRingCharacter(chi0);S:=CuspidalSubspace(ModularSymbols(chi,2,1));K:=BaseRing(S);d:=Dimension(S);dim+:=d*Degree(K);
 T:=HeckeOperator(S,qq);A:=(1+qq*Evaluate(chi,qq))*IdentityMatrix(K,d)-T;z:=Determinant(A);
 n:=Type(K) eq FldRat select Integers()!z else Integers()!Norm(z);NN*:=n;
 printf "NONIC247_FINITE_IMAGE_ODD q=%o order=%o dimension=%o coefficient_degree=%o determinant=%o norm=%o odd=%o time=%o\n",qq,Order(chi),d,Degree(K),Eltseq(z),n,IsOdd(n),Cputime(tt);
end for;
printf "NONIC247_FINITE_IMAGE_ODD_DONE q=%o dimension=%o odd=%o norm_product=%o time=%o\n",qq,dim,IsOdd(NN),NN,Cputime(tt);
if IsOdd(NN) then break;end if;
end for;quit;
