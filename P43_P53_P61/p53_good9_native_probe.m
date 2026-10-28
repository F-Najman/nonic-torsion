SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");
t0:=Cputime();p:=53;X:=MDX1(p,GF(2));C:=Curve(X);assert Genus(C) eq 92;
printf "NONIC53_MODEL genus=%o time=%o\n",Genus(C),Cputime(t0);
pl:=NoncuspidalPlaces(X,9);
printf "NONIC53_NC9 places=%o time=%o\n",#pl,Cputime(t0);
reps:=PlacesUpToDiamond(X,pl);
printf "NONIC53_NC9_REPS count=%o time=%o\n",#reps,Cputime(t0);
a:=PrimitiveRoot(p);
for i in [1..#reps] do
 P:=reps[i];D:=Divisor(P);T:=HeckeOperatorFast(X,3,P);
 A:=T-D-3*Divisor(DiamondOperator(X,3,P));
 B:=DiamondOperator(X,a,A)-A;assert Degree(B) eq 0;
 printf "NONIC53_TEST index=%o time=%o\n",i,Cputime(t0);
 flag:=IsPrincipal(B);
 printf "NONIC53_RESULT index=%o principal=%o time=%o\n",i,flag,Cputime(t0);
end for;
print "NONIC53_GOOD9_PROBE_COMPLETE";quit;
