SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
S:=ModularSymbolsH(89,[81],2,-1);M:=ModularSymbolsH(89,[81],2,0);MS:=CuspidalSubspace(M);e:=WindingElementProjection(MS);
A:=ChangeRing(IntegralHeckeOperator(S,3),Rationals());D:=(A^2-ChangeRing(IntegralHeckeOperator(S,9),Rationals()))/3;
Id:=Parent(A)!1;assert D^4 eq Id and D^2 ne Id;
p:=CharacteristicPolynomial(A);fac:=Factorization(p);assert [Degree(t[1]):t in fac] eq [1,1,5,6,14] and &and[t[2] eq 1:t in fac];
AM:=HeckeOperator(MS,3);flags:=[IsZero(e*Evaluate(ExactQuotient(p,f[1]),AM)):f in fac];assert flags eq [false,true,false,false,false];
mins:=[];for k in [1..#fac]do
 V:=Kernel(Evaluate(fac[k][1],A));B:=BasisMatrix(V);assert Nrows(B) eq Degree(fac[k][1]);R:=Solution(B,B*D);assert R*B eq B*D;
 mp:=MinimalPolynomial(R);Append(~mins,Eltseq(mp));if flags[k]then assert B*(D-Id) eq 0;end if;
 printf "NONIC89_WINDING_CHARACTER factor=%o dimension=%o T3_polynomial=%o winding_zero=%o diamond3_minpoly=%o time=%o\n",k,Nrows(B),fac[k][1],flags[k],mp,Cputime(tt);
end for;
assert Dimension(Kernel(D-Id)) eq 7 and Rank(D-Id) eq 20;
out:="P79_P89/data/P89/quotient27/finite_image_winding_complete_data.m";
PrintFile(out,Sprintf("N:=89;HGenerator:=81;Hecke3:=Matrix(Rationals(),27,27,%o);Diamond3:=Matrix(Rationals(),27,27,%o);MSHecke3:=Matrix(Rationals(),54,54,%o);MSWindingProjection:=%o;T3Factors:=%o;WindingZeroFlags:=%o;DiamondMinimalPolynomials:=%o;",Eltseq(A),Eltseq(D),Eltseq(AM),Eltseq(e),[Eltseq(f[1]):f in fac],flags,mins):Overwrite:=true);
PrintFile(out,"print \"NONIC89_FINITE_IMAGE_WINDING_DATA_READY\";");print "NONIC89_FINITE_IMAGE_ALPHA_CERTIFIED alpha=diamond3-1 rational_image_finite=true";quit;
