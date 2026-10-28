SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
S:=ModularSymbolsH(67,[38],2,-1);M:=ModularSymbolsH(67,[38],2,0);MS:=CuspidalSubspace(M);e:=WindingElementProjection(MS);
T3:=ChangeRing(IntegralHeckeOperator(S,3),Rationals());T7:=ChangeRing(IntegralHeckeOperator(S,7),Rationals());T2:=ChangeRing(IntegralHeckeOperator(S,2),Rationals());D:=(T2^2-ChangeRing(IntegralHeckeOperator(S,4),Rationals()))/2;Id:=Parent(D)!1;assert D^11 eq Id and D ne Id;
found:=false;
for c in [1..8]do
 A:=T3+c*T7;p:=CharacteristicPolynomial(A);fac:=Factorization(p);if exists{t:t in fac|t[2] ne 1}then continue;end if;
 AM:=HeckeOperator(MS,3)+c*HeckeOperator(MS,7);flags:=[IsZero(e*Evaluate(ExactQuotient(p,f[1]),AM)):f in fac];
 for power in [1]do
  alpha:=D^power-Id;ok:=true;for k in [1..#fac]do if flags[k]and BasisMatrix(Kernel(Evaluate(fac[k][1],A)))*alpha ne 0 then ok:=false;break;end if;end for;
  if ok then apow:=power;found:=true;break;end if;
 end for;if found then chosenC:=c;break;end if;
end for;assert found;mins:=[];
for k in [1..#fac]do
 V:=Kernel(Evaluate(fac[k][1],A));B:=BasisMatrix(V);assert Nrows(B) eq Degree(fac[k][1]);R:=Solution(B,B*D);assert R*B eq B*D;mp:=MinimalPolynomial(R);Append(~mins,Eltseq(mp));
 if flags[k]then assert B*(D^apow-Id) eq 0;end if;
 printf "NONIC67H3_WINDING_CHARACTER factor=%o dimension=%o splitting_polynomial=%o winding_zero=%o diamond2_minpoly=%o time=%o\n",k,Nrows(B),fac[k][1],flags[k],mp,Cputime(tt);
end for;
out:="P67/data/quotient45/finite_image_winding_data.m";
PrintFile(out,Sprintf("N:=67;HGenerator:=38;SplittingHecke7Coefficient:=%o;AlphaExponent:=%o;HeckeSplit:=Matrix(Rationals(),45,45,%o);Diamond2:=Matrix(Rationals(),45,45,%o);MSHeckeSplit:=Matrix(Rationals(),90,90,%o);MSWindingProjection:=%o;HeckeFactors:=%o;WindingZeroFlags:=%o;DiamondMinimalPolynomials:=%o;",chosenC,apow,Eltseq(A),Eltseq(D),Eltseq(AM),Eltseq(e),[Eltseq(f[1]):f in fac],flags,mins):Overwrite:=true);
PrintFile(out,"print \"NONIC67H3_FINITE_IMAGE_WINDING_DATA_READY\";");printf "NONIC67H3_FINITE_IMAGE_ALPHA_CERTIFIED sigma=diamond2 alpha=sigma^%o-1 rank=%o winding_zero_flags=%o splitting_operator=T3+%oT7 time=%o\n",apow,Rank(D^apow-Id),flags,chosenC,Cputime(tt);quit;
