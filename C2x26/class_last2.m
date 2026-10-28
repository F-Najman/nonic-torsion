/* Fresh proven genus25 rational-image query for ONLY the two exact
   degree9 source places surviving the complete66-orbit quotient filter. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);SetSeed(1);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
load "C2x26/data/N26_D9_threequotients_survivors.m";
assert #N26D9ThreeQuotientSurvivors eq 2;
X:=MDX11(2,26,GF(3));C:=Curve(X);assert Genus(C) eq 25;
cs:=Cusps(X);rat:=[c:c in cs|Degree(c) eq 1];assert #rat eq 18;c0:=rat[1];
fp:=FrobeniusPolynomial(JH(52,[27]),3);
print "C226_LAST2_PROVEN_CLASS_GROUP_START",Cputime(tt);
G,_,phi:=ClassGroupRationalReduction(C,c0:Proof:=true);
assert &*[a:a in Invariants(G)|a ne 0] eq Integers()!Evaluate(fp,1);
hg:=[phi(Divisor(c)-Divisor(c0)):c in rat];H:=sub<G|hg>;
assert Invariants(H) eq [2,14,266,3990,11970];
print "C226_LAST2_PROVEN_CLASS_GROUP",Invariants(G),Invariants(H),Cputime(tt);
R:=PolynomialRing(GF(3));records:=[* *];remaining:=[];
for rec in N26D9ThreeQuotientSurvivors do
 idx,pol,vals:=Explode(rec);K<z>:=ext<GF(3)|R!pol>;assert Degree(K) eq 9;
 bc:=[Evaluate(R!v,z):v in vals];x:=MDPlaceCommonZeros(X`_coordinates,bc);assert Degree(x) eq 9;
 E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);assert Order(Ls`P) eq 2 and Order(Ls`Q) eq 26 and Ls`P ne 13*Ls`Q;
 cl:=phi(Divisor(x)-9*Divisor(c0));inside:=cl in H;Append(~records,[*idx,pol,vals,Eltseq(cl),inside*]);
 if inside then Append(~remaining,idx);end if;
 print "C226_LAST2_PROVEN_CLASS_CASE",idx,Eltseq(cl),inside,Cputime(tt);
 Write("C2x26/data/class_last2_partial_arithmetic.txt",Sprint([*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(g):g in hg],records,remaining,Coefficients(fp)*]));
end for;
Write("C2x26/data/class_last2_arithmetic.txt",Sprint([*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(g):g in hg],records,remaining,Coefficients(fp)*]):Overwrite:=true);
printf "C226_LAST2_PROVEN_CLASS_DONE tested=2 remaining=%o complete=true cpu=%o\n",#remaining,Cputime(tt);quit;
