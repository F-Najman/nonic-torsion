/* Pure degree9 supports surviving all three X1(26) rational-image maps.
   Every source automorphism commutes with odd Hecke correspondences. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
load "C2x26/data/N26_D9_threequotients_survivors.m";
assert #N26D9ThreeQuotientSurvivors eq 2;
X:=MDX11(2,26,GF(3));assert Genus(X) eq 25;R:=PolynomialRing(GF(3));records:=[* *];remaining:=[];
for rec in N26D9ThreeQuotientSurvivors do
 idx,pol,vals:=Explode(rec);K<z>:=ext<GF(3)|R!pol>;assert Degree(K) eq 9;
 bc:=[Evaluate(R!v,z):v in vals];x:=MDPlaceCommonZeros(X`_coordinates,bc);assert Degree(x) eq 9;
 killed:=false;
 for q in [5,7,11] do
  T:=HeckeOperatorFast(X,q,x);delta:=DiamondOperator(X,q,x);A:=T-Divisor(x)-q*Divisor(delta);
  assert Degree(T) eq 9*(q+1) and Degree(A) eq 0;
  principal:=IsPrincipal(A);Append(~records,[*idx,q,principal*]);
  print "N26_D9_LAST2_HECKE",idx,q,principal,Cputime(tt);
  if q eq 5 and idx eq N26D9ThreeQuotientSurvivors[1][1] then
   E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);ys:=[];ms:=[];
   for ph in MDIsogenies(E,5) do
    LL:=rec<recformat<P:PtEll,Q:PtEll>|P:=ph(Ls`P),Q:=ph(Ls`Q)>;
    yy:=ModuliPointFast(X,Codomain(ph),LL);i:=Index(ys,yy);
    if i eq 0 then Append(~ys,yy);Append(~ms,1);else ms[i]+:=1;end if;
   end for;
   reference:=DivisorGroup(Curve(X))!0;
   for i in [1..#ys] do mult:=ms[i]*9/Degree(ys[i]);assert Denominator(mult) eq 1;reference+:=Integers()!mult*Divisor(ys[i]);end for;
   assert reference eq T;print "N26_D9_LAST2_FULL_NATIVE_T5_EQUAL",Degree(T),Cputime(tt);
  end if;
  if not principal then killed:=true;break;end if;
 end for;
 if not killed then Append(~remaining,idx);end if;
end for;
Write("C2x26/data/N26_D9_last2_Hecke_arithmetic.txt",Sprint([*records,remaining*]):Overwrite:=true);
printf "N26_D9_LAST2_HECKE_DONE tested=2 remaining=%o complete=true time=%o\n",#remaining,Cputime(tt);quit;
