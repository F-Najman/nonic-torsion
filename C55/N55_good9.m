SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");Attach("common/fast_hecke_cyclic_commonzero.m");
X:=MDX1(55,GF(2));P9:=NoncuspidalPlaces(X,9);rs:=PlacesUpToDiamondFast(X,P9);
printf "NONIC55_GOOD9_SETUP genus=%o all_places=%o representatives=%o time=%o\n",Genus(X),#P9,#rs,Cputime(t0);
for k in [1..#rs] do
 D:=Divisor(rs[k]);T:=HeckeOperatorCommonZeros(X,3,D);A:=T-D-3*DiamondOperator(X,3,D);assert Degree(A) eq 0;
 flag:=IsPrincipal(A);printf "NONIC55_GOOD9_CASE index=%o principal=%o time=%o\n",k,flag,Cputime(t0);assert not flag;
 if k eq 1 then Native:=HeckeOperator(X,3,D);assert Native eq T;printf "NONIC55_FULL_NATIVE_EQUAL degree=%o time=%o\n",Degree(T),Cputime(t0);end if;
end for;
printf "NONIC55_GOOD9_DONE all_places=%o representatives=%o complete=true time=%o\n",#P9,#rs,Cputime(t0);quit;
