SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
load "C2x28/data/four_quotients_survivors.m";
rr:=[r:r in C228FourQuotientSurvivors|r[1] eq 9 and r[2] eq 10][1];d,i,pol,vals:=Explode(rr);R:=PolynomialRing(GF(3));K<z>:=ext<GF(3)|R!pol>;bc:=[Evaluate(R!v,z):v in vals];
X:=MDX11(2,28,GF(3));x:=MDPlaceCommonZeros(X`_coordinates,bc);assert Degree(x) eq 9;
T:=HeckeOperatorFast(X,5,x);E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);ys:=[];ms:=[];
for ph in MDIsogenies(E,5) do
 LL:=rec<recformat<P:PtEll,Q:PtEll>|P:=ph(Ls`P),Q:=ph(Ls`Q)>;yy:=ModuliPointFast(X,Codomain(ph),LL);j:=Index(ys,yy);
 if j eq 0 then Append(~ys,yy);Append(~ms,1);else ms[j]+:=1;end if;
end for;
reference:=DivisorGroup(Curve(X))!0;
for j in [1..#ys] do m:=ms[j]*9/Degree(ys[j]);assert Denominator(m) eq 1;reference+:=Integers()!m*Divisor(ys[j]);end for;
assert reference eq T and Degree(T) eq 54;
A:=T-Divisor(x)-5*Divisor(DiamondOperator(X,5,x));assert Degree(A) eq 0;pr:=IsPrincipal(A);assert not pr;
print "C228_FULL_NATIVE_T5_CONTROL_EQUAL degree=54 nonprincipal=true complete=true",Cputime(tt);quit;
