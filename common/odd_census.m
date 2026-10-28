SetColumns(0);SetQuitOnError(true);SetMemoryLimit(10^9);tt:=Cputime();
if Type(N)eq MonStgElt then N:=StringToInteger(N);end if;
if Type(FieldDegree)eq MonStgElt then FieldDegree:=StringToInteger(FieldDegree);end if;
assert N in[57,65,85,91,93,117] and FieldDegree in[6..9];F<a>:=GF(2^FieldDegree);
Code:=func<v|[Integers()!c:c in Eltseq(v,GF(2))]>;
out:=Open(OutFile,"w");fprintf out,"N:=%o;FieldDegree:=%o;FieldPolynomial:=%o;Records:=[*\n",N,FieldDegree,[Integers()!c:c in Eltseq(DefiningPolynomial(F,GF(2)))];
nt:=0;divisible:=0;qualified:=0;js:={};jr:={};first:=true;
for j in F do
 for E in Twists(EllipticCurveFromjInvariant(j))do
  nt+:=1;order:=#E;if order mod N ne 0 then continue;end if;divisible+:=1;
  A,mp:=AbelianGroup(E);iv:=Invariants(A);ex:=iv[#iv];if ex mod N ne 0 then continue;end if;
  P:=(ex div N)*mp(A.#iv);assert N*P eq E!0 and &and[(N div p)*P ne E!0:p in PrimeDivisors(N)];qualified+:=1;Include(~js,j);Include(~jr,MinimalPolynomial(j,GF(2)));
  if not first then fprintf out,",\n";end if;first:=false;
  fprintf out,"<%o,%o,%o,%o,%o,%o>",Code(j),[Code(v):v in aInvariants(E)],[Code(P[1]),Code(P[2])],order,iv,[Integers()!c:c in Eltseq(MinimalPolynomial(j,GF(2)))];
 end for;
end for;
fprintf out,"\n*];TwistsChecked:=%o;OrderDivisibleCurves:=%o;QualifiedCurves:=%o;QualifiedJValues:=%o;QualifiedJPolynomials:=%o;Complete:=true;\n",nt,divisible,qualified,#js,[[Integers()!c:c in Eltseq(p)]:p in jr];delete out;
printf "NONIC_ODD_CENSUS_COMPLETE N=%o field_degree=%o twists=%o divisible=%o qualified=%o j_values=%o j_orbit_degrees=%o time=%o\n",N,FieldDegree,nt,divisible,qualified,#js,Sort([Degree(p):p in jr]),Cputime(tt);quit;
