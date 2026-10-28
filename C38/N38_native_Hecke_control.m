SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");Attach("common/fast_hecke_cyclic_commonzero.m");
X:=MDX1(38,GF(3));S:=NoncuspidalPlaces(X,9);assert #S eq 2127;P:=S[1];D:=Divisor(P);
printf "NONIC38_NATIVE_CONTROL_INPUT degree=%o totalNC9=%o time=%o\n",Degree(P),#S,Cputime(t0);
q:=5;F:=HeckeOperatorCommonZeros(X,q,D);T:=HeckeOperator(X,q,D);
assert F eq T;B:=T-D-q*DiamondOperator(X,q,D);flag:=IsPrincipal(B);
printf "NONIC38_NATIVE_HECKE_EQUAL q=5 degree=%o principal=%o time=%o\n",Degree(T),flag,Cputime(t0);
print "NONIC38_NATIVE_CONTROL_COMPLETE";quit;
