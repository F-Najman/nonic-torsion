SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();N:=65;HH:=[u:u in[1..64]|GCD(u,65)eq 1 and KroneckerSymbol(65,u)eq 1];assert #HH eq 24;
S:=CuspidalSubspace(ModularSymbolsH(N,HH,2,1));g:=Dimension(S);T3:=HeckeOperator(S,3);D3:=(T3^2-HeckeOperator(S,9))/3;I:=IdentityMatrix(Rationals(),g);assert D3^2 eq I;
printf "NONIC65_QUOTIENT_PREFLIGHT genus=%o minus_rank=%o time=%o\n",g,Rank(D3-I),Cputime(tt);decomp:=NewformDecomposition(S);minus:=[];
for A in decomp do ff:=qEigenform(A,30);cc:=(Coefficient(ff,3)^2-Coefficient(ff,9))/3;printf "NONIC65_QUOTIENT_FACTOR dim=%o level=%o delta3=%o winding_nonzero=%o\n",Dimension(A),Level(AssociatedNewSpace(A)),cc,not IsZero(WindingElementProjection(A));if cc eq -1 then Append(~minus,A);end if;end for;
orders:=[];
for ell in[2,3,7,11,17,19,23,29,31]do
 v:=1;
 for A in minus do n:=Dimension(A);AA:=HeckeOperator(A,ell);DD:=(AA^2-HeckeOperator(A,ell^2))/ell;v*:=Integers()!Determinant(IdentityMatrix(Rationals(),n)-AA+ell*DD);end for;
 Append(~orders,<ell,v>);printf "NONIC65_PRYM_ORDER prime=%o order=%o gcd=%o time=%o\n",ell,v,GCD([v[2]:v in orders]),Cputime(tt);
end for;
Write("C65/data/quotient9/prym_preflight.txt",Sprint([*g,Rank(D3-I),orders,GCD([v[2]:v in orders]) *]):Overwrite:=true);print "NONIC65_QUOTIENT_PREFLIGHT_COMPLETE";quit;
