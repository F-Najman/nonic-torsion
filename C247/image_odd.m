SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();NN:=1;dim:=0;
for chi0 in GaloisConjugacyRepresentatives(FullDirichletGroup(247)) do
 if IsOdd(chi0) or Evaluate(chi0,40) ne 1 or Order(chi0) eq 1 then continue;end if;
 chi:=MinimalBaseRingCharacter(chi0);S:=CuspidalSubspace(ModularSymbols(chi,2,1));K:=BaseRing(S);d:=Dimension(S);dim+:=d*Degree(K);
 T:=HeckeOperator(S,3);A:=(1+3*Evaluate(chi,3))*IdentityMatrix(K,d)-T;z:=Determinant(A);
 n:=Type(K) eq FldRat select Integers()!z else Integers()!Norm(z);NN*:=n;
 printf "NONIC247_FINITE_IMAGE_F3 order=%o dimension=%o coefficient_degree=%o determinant=%o norm=%o odd=%o time=%o\n",Order(chi),d,Degree(K),Eltseq(z),n,IsOdd(n),Cputime(tt);
end for;
printf "NONIC247_FINITE_IMAGE_F3_DONE dimension=%o odd=%o norm_product=%o time=%o\n",dim,IsOdd(NN),NN,Cputime(tt);quit;
