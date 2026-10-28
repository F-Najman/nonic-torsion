SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
S:=ModularSymbolsH(79,[27],2,-1);M:=ModularSymbolsH(79,[27],2,0);MS:=CuspidalSubspace(M);e:=WindingElementProjection(MS);
A:=ChangeRing(IntegralHeckeOperator(S,3),Rationals());D:=(A^2-ChangeRing(IntegralHeckeOperator(S,9),Rationals()))/3;
Id:=Parent(A)!1;assert D^3 eq Id and D ne Id;
p:=CharacteristicPolynomial(A);fac:=Factorization(p);assert [Degree(t[1]):t in fac] eq [1,5,12] and &and[t[2] eq 1:t in fac];
AM:=HeckeOperator(MS,3);flags:=[IsZero(e*Evaluate(ExactQuotient(p,f[1]),AM)):f in fac];assert flags eq [true,false,false];
mins:=[];for k in [1..#fac]do
 V:=Kernel(Evaluate(fac[k][1],A));B:=BasisMatrix(V);assert Nrows(B) eq Degree(fac[k][1]);R:=Solution(B,B*D);assert R*B eq B*D;
 mp:=MinimalPolynomial(R);Append(~mins,Eltseq(mp));if flags[k]then assert B*(D-Id) eq 0;end if;
 printf "NONIC79_WINDING_CHARACTER factor=%o dimension=%o T3_polynomial=%o winding_zero=%o diamond3_minpoly=%o time=%o\n",k,Nrows(B),fac[k][1],flags[k],mp,Cputime(tt);
end for;
assert Dimension(Kernel(D-Id)) eq 6 and Rank(D-Id) eq 12;
out:="P79_P89/data/P79/quotient18/finite_image_winding_complete_data.m";
PrintFile(out,Sprintf("N:=79;HGenerator:=27;Hecke3:=Matrix(Rationals(),18,18,%o);Diamond3:=Matrix(Rationals(),18,18,%o);MSHecke3:=Matrix(Rationals(),36,36,%o);MSWindingProjection:=%o;T3Factors:=%o;WindingZeroFlags:=%o;DiamondMinimalPolynomials:=%o;",Eltseq(A),Eltseq(D),Eltseq(AM),Eltseq(e),[Eltseq(f[1]):f in fac],flags,mins):Overwrite:=true);
PrintFile(out,"print \"NONIC79_FINITE_IMAGE_WINDING_DATA_READY\";");print "NONIC79_FINITE_IMAGE_ALPHA_CERTIFIED alpha=diamond3-1 rational_image_finite=true";quit;
