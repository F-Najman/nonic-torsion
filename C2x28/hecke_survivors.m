/* Degree9 surviving single supports; all cubic cusp corrections for d6.
   d>=7 has cuspcomplement<=2 and every such closed cusp is Q-defined. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
load "C2x28/data/four_quotients_survivors.m";
assert #C228FourQuotientSurvivors eq 10;
X:=MDX11(2,28,GF(3));C:=Curve(X);FF:=FunctionField(C);cs:=Cusps(X);cub:=[c:c in cs|Degree(c) eq 3];assert #cub eq 4;
EE:=EllipticCurve([FF!f:f in X`_E]);QQ:=EE![FF!f:f in X`_Q];bb:=QQ[1];cc:=aInvariants(EE)[2];dd:=aInvariants(EE)[4];assert dd eq bb*(1-bb-cc);
fs:=[bb,dd,cc^2-4*dd];for k in [1..14] do Append(~fs,(k*QQ)[1]);end for;
sigs:=[[Valuation(f,c):f in fs]:c in cub];us:=Sort(Setseq(Seqset(sigs)));assert #us eq 2;
pairs:=[[i:i in [1..4]|sigs[i] eq s]:s in us];assert &and[#p eq 2:p in pairs];vs:=[];
for p in pairs do for eps in [-1,1] do v:=[0:i in [1..4]];v[p[1]]:=-eps;v[p[2]]:=eps;Append(~vs,v);end for;end for;
R:=PolynomialRing(GF(3));tests:=[* *];remaining:=[* *];
for rec in C228FourQuotientSurvivors do
 d,i,pol,vals:=Explode(rec);assert d in [6..9];K<z>:=ext<GF(3)|R!pol>;bc:=[Evaluate(R!v,z):v in vals];x:=MDPlaceCommonZeros(X`_coordinates,bc);assert Degree(x) eq d;killed:=false;
 for q in [5,11] do
  T:=HeckeOperatorFast(X,q,x);A:=T-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));assert Degree(T) eq d*(q+1) and Degree(A) eq 0;
  corrections:=d eq 6 select [[0,0,0,0]] cat [[v[k]+q*w[k]:k in [1..4]]:v,w in vs] else [[0,0,0,0]];
  assert #Seqset(corrections) eq (d eq 6 select 17 else 1);flags:=[];
  for v in corrections do
   B:=A+&+[v[k]*Divisor(cub[k]):k in [1..4]];assert Degree(B) eq 0;Append(~flags,IsPrincipal(B) select 1 else 0);
  end for;
  Append(~tests,[*d,i,q,corrections,flags*]);print "C228_D9_HECKE_CASE",d,i,q,&+flags,Cputime(tt);
  if &+flags eq 0 then killed:=true;break;end if;
 end for;
 if not killed then Append(~remaining,rec);end if;
 Write("C2x28/data/hecke_survivors_arithmetic.txt",Sprint([*sigs,pairs,vs,tests,remaining*]):Overwrite:=true);
 Write("C2x28/data/hecke_survivors_remaining.m","C228HeckeSurvivors:=" cat Sprint(remaining) cat ";":Overwrite:=true);
end for;
printf "C228_D9_HECKE_DONE tested=10 remaining=%o complete=true cpu=%o\n",#remaining,Cputime(tt);quit;
