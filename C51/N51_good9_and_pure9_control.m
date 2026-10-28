SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");Attach("common/fast_hecke_cyclic_commonzero.m");
X:=MDX1(51,GF(2));C:=Curve(X);FF:=FunctionField(C);AF,afmap:=AlgorithmicFunctionField(FF);g:=Genus(AF);assert g eq 65;
printf "NONIC51_MODEL genus=65 time=%o\n",Cputime(t0);
O,om:=SpaceOfHolomorphicDifferentials(AF);B:=[om(v):v in Basis(O)];CM:=Matrix(GF(2),[Eltseq(Cartier(w) @@ om):w in B]);U:=RowSpace(CM-IdentityMatrix(GF(2),g));assert Dimension(U) eq 62;
jaf:=afmap(jInvariant(EllipticCurve(X`_E)));cs:=Poles(jaf);assert {*Degree(c):c in cs*} eq {*1^^16,2^^8,8^^4*};rat:=[c:c in cs|Degree(c) eq 1];
printf "NONIC51_CUSPS_DONE rational=%o time=%o\n",#rat,Cputime(t0);
V,vm:=DifferentialSpace(9*Divisor(rat[1]));W:=sub<VectorSpace(GF(2),g)|[Eltseq(vm(v) @@ om):v in Basis(V)]>;rank9:=Dimension(U)-Dimension(U meet W);assert rank9 eq 9;
printf "NONIC51_PURE9_CUSP_CONTROL rank=%o covered=16 dimensionU=%o time=%o\n",rank9,Dimension(U),Cputime(t0);
P9:=NoncuspidalPlaces(X,9);rs:=PlacesUpToDiamondFast(X,P9);
printf "NONIC51_GOOD9_SETUP all_places=%o representatives=%o time=%o\n",#P9,#rs,Cputime(t0);
for k in [1..#rs] do
 P:=rs[k];D:=Divisor(P);killed:=false;
 for q in [5,7,11,13,19] do
  T:=HeckeOperatorCommonZeros(X,q,D);A:=T-D-q*DiamondOperator(X,q,D);assert Degree(A) eq 0;
  flag:=IsPrincipal(A);printf "NONIC51_GOOD9_CASE index=%o q=%o principal=%o time=%o\n",k,q,flag,Cputime(t0);
  if not flag then killed:=true;break;end if;
 end for;
 assert killed;
end for;
printf "NONIC51_GOOD9_AND_PURE9_DONE all_places=%o representatives=%o complete=true time=%o\n",#P9,#rs,Cputime(t0);quit;
