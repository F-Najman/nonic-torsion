SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
N:=65;h:=0;gg:=9;pp:=2;prec:=1100;small:=120;
HH:=[u:u in[1..64]|GCD(u,65)eq 1 and KroneckerSymbol(65,u)eq 1];assert #HH eq 24;
MS:=CuspidalSubspace(ModularSymbolsH(N,HH,2,1));assert Dimension(MS) eq gg;
Dcmp:=NewformDecomposition(MS);
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
rows:=[];constits:=[];
for A in Dcmp do
 lev:=Level(AssociatedNewSpace(A));assert N mod lev eq 0;
 ff:=qEigenform(A,prec);KF:=BaseRing(Parent(ff));block:=[];
 if KF cmpeq Rationals() then Append(~block,[Rationals()!Coefficient(ff,n):n in [0..prec-1]]);else
  dd:=#Flat(KF!1);co:=[Flat(KF!Coefficient(ff,n)):n in [0..prec-1]];assert &and[#v eq dd:v in co];
  for k in [1..dd] do Append(~block,[v[k]:v in co]);end for;
 end if;
 for d in Divisors(N div lev) do
  for r in block do Append(~rows,[n mod d eq 0 select r[n div d+1] else Rationals()!0:n in [0..prec-1]]);end for;
 end for;
 Append(~constits,<lev,#block,Divisors(N div lev),Dimension(A)>);
 printf "NONIC65_CONSTITUENT level=%o coefficient_dimension=%o old_divisors=%o symbol_dimension=%o time=%o\n",lev,#block,Divisors(N div lev),Dimension(A),Cputime(tt);
end for;assert #rows eq gg;
Qraw:=Matrix(Rationals(),rows);Q0:=Submatrix(Qraw,1,1,gg,small);assert Rank(Q0) eq gg;
den:=LCM([Denominator(a):a in Eltseq(Q0)]);QS:=ChangeRing(Saturation(ChangeRing(den*Q0,Integers())),Rationals());assert Nrows(QS) eq gg and RowSpace(QS) eq RowSpace(Q0);
Tr:=Solution(Q0,QS);full:=Tr*Qraw;assert Submatrix(full,1,1,gg,small) eq QS;
assert &and[Denominator(a) mod pp ne 0:a in Eltseq(full)];F5:=GF(pp);Q5:=ChangeRing(full,F5);assert Rank(Submatrix(Q5,1,1,gg,small)) eq gg;
E,U:=EchelonForm(Submatrix(Q5,1,1,gg,small));assert U*Submatrix(Q5,1,1,gg,small) eq E;
Q5:=U*Q5;
PS<q>:=PowerSeriesRing(F5,prec);Q:=[PS!Eltseq(v)+O(q^prec):v in Rows(Q5)];
pairs:=[];for i in [1..gg],j in [i..gg] do Append(~pairs,<i,j>);end for;
A:=Matrix(F5,[[Coefficient(Q[z[1]]*Q[z[2]],n):n in [0..small-1]]:z in pairs]);assert Rank(A) eq 3*gg-3;rels:=BasisMatrix(Kernel(A));assert Nrows(rels) eq 21;
out:="C65/data/quotient9/canonical_data.m";
PrintFile(out,Sprintf("N:=65;genus:=9;p:=2;H_generator:=0;precision:=%o;SmallPrecision:=%o;Constituents:=%o;qcoeffs:=Matrix(GF(2),9,%o,%o);pairs:=%o;QuadRelations:=Matrix(GF(2),21,45,%o);",prec,small,constits,prec,Eltseq(Q5),pairs,Eltseq(rels)):Overwrite:=true);
P:=ProjectiveSpace(F5,gg-1);R:=CoordinateRing(P);mons:=[R.(z[1])*R.(z[2]):z in pairs];eqs:=[&+[v[i]*mons[i]:i in [1..#pairs]]:v in Rows(rels)];C:=Curve(P,eqs);assert Dimension(C) eq 1 and Degree(C) eq 16 and ArithmeticGenus(C) eq gg;
printf "NONIC65_CANONICAL_MODEL_READY genus9=true degree16=true saturated_mod2=true all_old_divisors=true time=%o\n",Cputime(tt);
LS<qq>:=LaurentSeriesRing(F5,prec);F:=[LS!Eltseq(v)+O(qq^prec):v in Rows(Q5)];
foundTriple:=false;
for i in [gg..1 by -1] do
 for j in [gg..1 by -1] do
  if i eq j then continue;end if;
  if not IsZero(Derivative(F[j]/F[i])) then
   k:=[k:k in [gg..1 by -1]|k ne i and k ne j][1];ix:=[i,j,k];foundTriple:=true;break;
  end if;
 end for;
 if foundTriple then break;end if;
end for;assert foundTriple;
a,bb,cc:=Explode([F[i]:i in ix]);uu:=bb/a;vv:=cc/a;assert not IsZero(Derivative(uu));
Dmax:=2*gg-1-Min([Valuation(w):w in [a,bb,cc]]);assert 18*Dmax+1 lt prec;
printf "NONIC65_PLANE_PROJECTION orders=%o maxdegree=%o separable=true time=%o\n",[Valuation(w):w in [a,bb,cc]],Dmax,Cputime(tt);
PP<x,y,z>:=ProjectiveSpace(F5,2);RR:=CoordinateRing(PP);found:=false;
for dd in [2..Dmax] do
 mm:=Setseq(MonomialsOfDegree(RR,dd));ser:=[Evaluate(m,[bb,cc,a]):m in mm];AA:=Matrix(F5,[[Coefficient(s,n):n in [0..18*dd]]:s in ser]);KK:=Kernel(AA);
 if Dimension(KK) gt 0 then assert Dimension(KK) eq 1;v:=KK.1;peq:=&+[v[i]*mm[i]:i in [1..#mm]];assert IsIrreducible(peq);planeDeg:=dd;found:=true;break;end if;
end for;assert found;
PrintFile("C65/data/quotient9/plane_equation.m",Sprintf("N:=65;genus:=9;p:=2;plane_degree:=%o;PlaneExponents:=%o;PlaneCoefficients:=%o;ProjectionIndices:=%o;",planeDeg,[Exponents(m):m in Monomials(peq)],[Integers()!v:v in Coefficients(peq)],ix):Overwrite:=true);
printf "NONIC65_PLANE_EQUATION_READY degree=%o time=%o\n",planeDeg,Cputime(tt);
AA<xx,yy>:=AffineSpace(F5,2);RR2:=CoordinateRing(AA);CA:=Curve(AA,Evaluate(peq,[RR2.1,RR2.2,RR2!1]));CC:=ProjectiveClosure(CA);FF:=FunctionField(CC);AF,am:=AlgorithmicFunctionField(FF);gx:=Genus(AF);assert gx eq gg;
uAF:=am(FF.1);vAF:=am(FF.2);KAF:=BaseField(AF);
if uAF eq AF!KAF.1 and vAF eq AF.1 then baseSeries:=uu;genSeries:=vv;baseCoord:=1;genCoord:=2;
elif vAF eq AF!KAF.1 and uAF eq AF.1 then baseSeries:=vv;genSeries:=uu;baseCoord:=2;genCoord:=1;
else error "Unexpected AF coordinate change";end if;
printf "NONIC65_PLANE_NORMALIZATION_READY genus=%o AFdegree=%o time=%o\n",gx,Degree(AF),Cputime(tt);
Om,om:=SpaceOfHolomorphicDifferentials(AF);B:=[om(w):w in Basis(Om)];eta:=Differential(uAF);BC:=[w/eta:w in B];
function EvalAF(w)
 co:=Eltseq(w);return &+[Evaluate(Numerator(co[j+1]),baseSeries)/Evaluate(Denominator(co[j+1]),baseSeries)*genSeries^j:j in [0..#co-1]];
end function;
forms:=[qq*EvalAF(w)*Derivative(uu):w in BC];assert &and[AbsolutePrecision(f) ge small and Valuation(f) ge 1:f in forms];OM:=Matrix(F5,[[Coefficient(f,n):n in [0..small-1]]:f in forms]);assert Rank(OM) eq gg;TC:=Solution(OM,Submatrix(Q5,1,1,gg,small));assert Rank(TC) eq gg;
Z:=[&+[TC[i,j]*BC[j]:j in [1..gg]]:i in [1..gg]];cz:=[Z[p[1]]*Z[p[2]]:p in pairs];assert &and[&+[v[i]*cz[i]:i in [1..#pairs]] eq 0:v in Rows(rels)];
la:=Z[ix[1]];assert Z[ix[2]]/la eq uAF and Z[ix[3]]/la eq vAF;
function CodeK(w)return <[Integers()!a:a in Eltseq(Numerator(w))],[Integers()!a:a in Eltseq(Denominator(w))]>;end function;
function CodeAF(w)return [CodeK(c):c in Eltseq(w)];end function;
PrintFile("C65/data/quotient9/plane_model_data.m",Sprintf("N:=65;genus:=9;p:=2;AF_degree:=%o;AFBaseCoordinate:=%o;AFGeneratorCoordinate:=%o;PlaneExponents:=%o;PlaneCoefficients:=%o;ProjectionIndices:=%o;AFPolynomialData:=%o;CanonicalCoordinateData:=%o;",Degree(AF),baseCoord,genCoord,[Exponents(m):m in Monomials(peq)],[Integers()!v:v in Coefficients(peq)],ix,[CodeK(w):w in Eltseq(DefiningPolynomial(AF))],[CodeAF(w):w in Z]):Overwrite:=true);
printf "NONIC65_MODEL_CERTIFIED genus9=true all21canonical_quadrics=true inverse_projection=true time=%o\n",Cputime(tt);quit;
