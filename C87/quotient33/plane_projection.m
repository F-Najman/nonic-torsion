SetColumns(0);SetQuitOnError(true);SetMemoryLimit(32*10^9);tt:=Cputime();
load "C87/data/quotient33/canonical_data.m";
gg:=33;prec:=3000;MS:=CuspidalSubspace(ModularSymbolsH(87,[16],2,1));decomp:=NewformDecomposition(MS);
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
rows:=[];
for A in decomp do
 lev:=Level(AssociatedNewSpace(A));assert lev in [29,87];ff:=qEigenform(A,prec);KF:=BaseRing(Parent(ff));first:=#rows+1;
 if KF cmpeq Rationals() then Append(~rows,[Rationals()!Coefficient(ff,n):n in [0..prec-1]]);else
  dd:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(ff,n)):n in [0..prec-1]];
  for k in [1..dd] do Append(~rows,[z[k]:z in cofs]);end for;
 end if;
 if lev eq 29 then oldRows:=rows[first..#rows];for r in oldRows do Append(~rows,[n mod 3 eq 0 select r[n div 3+1] else Rationals()!0:n in [0..prec-1]]);end for;end if;
end for;assert #rows eq gg;
Q0:=Matrix(Rationals(),rows);old:=ChangeRing(qcoeffs,Rationals());T:=Solution(Submatrix(Q0,1,1,gg,190),old);full:=T*Q0;assert Submatrix(full,1,1,gg,190) eq old;
assert &and[IsOdd(Denominator(a)):a in Eltseq(full)];Q2:=ChangeRing(full,GF(2));E,U:=EchelonForm(Submatrix(Q2,1,1,gg,190));assert U*Submatrix(Q2,1,1,gg,190) eq E;
LS<q>:=LaurentSeriesRing(GF(2),prec);F:=[LS!Eltseq(v)+O(q^prec):v in Rows(U*Q2)];ix:=[31,32,33];a,bb,cc:=Explode([F[i]:i in ix]);uu:=bb/a;vv:=cc/a;
if IsZero(Derivative(uu)) then assert not IsZero(Derivative(vv));ix:=[31,33,32];a,bb,cc:=Explode([F[i]:i in ix]);uu:=bb/a;vv:=cc/a;end if;assert not IsZero(Derivative(uu));
Dmax:=65-Valuation(a);assert 80*Dmax+1 lt prec;
printf "NONIC87_PLANE_PROJECTION orders=%o maxdegree=%o separable=true time=%o\n",[Valuation(w):w in [a,bb,cc]],Dmax,Cputime(tt);
P<x,y,z>:=ProjectiveSpace(GF(2),2);R:=CoordinateRing(P);found:=false;
for dd in [2..Dmax] do
 mons:=Setseq(MonomialsOfDegree(R,dd));ser:=[Evaluate(m,[bb,cc,a]):m in mons];A:=Matrix(GF(2),[[Coefficient(s,n):n in [0..80*dd]]:s in ser]);K:=Kernel(A);
 if Dimension(K) gt 0 then assert Dimension(K) eq 1;v:=K.1;peq:=&+[v[i]*mons[i]:i in [1..#mons]];assert IsIrreducible(peq);planeDeg:=dd;found:=true;break;end if;
end for;assert found;
proj:=Matrix(GF(2),[Eltseq(U[i]):i in ix]);PrintFile("C87/data/quotient33/plane_equation.m",Sprintf("N:=87;plane_degree:=%o;PlaneExponents:=%o;ProjectionRows:=Matrix(GF(2),3,33,%o);",planeDeg,[Exponents(m):m in Monomials(peq)],Eltseq(proj)):Overwrite:=true);
printf "NONIC87_PLANE_EQUATION_READY degree=%o time=%o\n",planeDeg,Cputime(tt);
AA<xx,yy>:=AffineSpace(GF(2),2);R2:=CoordinateRing(AA);CA:=Curve(AA,Evaluate(peq,[R2.1,R2.2,R2!1]));C:=ProjectiveClosure(CA);FF:=FunctionField(C);AF,am:=AlgorithmicFunctionField(FF);gx:=Genus(AF);assert gx eq gg;
uAF:=am(FF.1);vAF:=am(FF.2);KAF:=BaseField(AF);
if uAF eq AF!KAF.1 and vAF eq AF.1 then baseSeries:=uu;genSeries:=vv;baseCoord:=1;genCoord:=2;
elif vAF eq AF!KAF.1 and uAF eq AF.1 then baseSeries:=vv;genSeries:=uu;baseCoord:=2;genCoord:=1;
else error "Unexpected AF coordinate change";end if;
printf "NONIC87_PLANE_NORMALIZATION_READY genus=%o AFdegree=%o time=%o\n",gx,Degree(AF),Cputime(tt);
Om,om:=SpaceOfHolomorphicDifferentials(AF);B:=[om(w):w in Basis(Om)];eta:=Differential(uAF);BC:=[w/eta:w in B];
function EvalAF(w)
 co:=Eltseq(w);return &+[Evaluate(Numerator(co[j+1]),baseSeries)/Evaluate(Denominator(co[j+1]),baseSeries)*genSeries^j:j in [0..#co-1]];
end function;
forms:=[q*EvalAF(w)*Derivative(uu):w in BC];assert &and[AbsolutePrecision(f) ge 190 and Valuation(f) ge 1:f in forms];OM:=Matrix(GF(2),[[Coefficient(f,n):n in [0..189]]:f in forms]);assert Rank(OM) eq gg;TC:=Solution(OM,ChangeRing(qcoeffs,GF(2)));assert Rank(TC) eq gg;
Z:=[&+[TC[i,j]*BC[j]:j in [1..gg]]:i in [1..gg]];cz:=[Z[p[1]]*Z[p[2]]:p in pairs];assert &and[&+[v[i]*cz[i]:i in [1..#pairs]] eq 0:v in Rows(QuadRelations)];
la:=&+[proj[1,i]*Z[i]:i in [1..gg]];assert (&+[proj[2,i]*Z[i]:i in [1..gg]])/la eq uAF and (&+[proj[3,i]*Z[i]:i in [1..gg]])/la eq vAF;
function CodeK(w)return <[Integers()!a:a in Eltseq(Numerator(w))],[Integers()!a:a in Eltseq(Denominator(w))]>;end function;
function CodeAF(w)return [CodeK(c):c in Eltseq(w)];end function;
PrintFile("C87/data/quotient33/plane_model_data.m",Sprintf("N:=87;genus:=33;AF_degree:=%o;AFBaseCoordinate:=%o;AFGeneratorCoordinate:=%o;PlaneExponents:=%o;ProjectionRows:=Matrix(GF(2),3,33,%o);AFPolynomialData:=%o;CanonicalCoordinateData:=%o;",Degree(AF),baseCoord,genCoord,[Exponents(m):m in Monomials(peq)],Eltseq(proj),[CodeK(w):w in Eltseq(DefiningPolynomial(AF))],[CodeAF(w):w in Z]):Overwrite:=true);
printf "NONIC87_QUOTIENT_PLANE_CERTIFIED genus33=true all465canonical_quadrics=true inverse_projection=true time=%o\n",Cputime(tt);quit;
