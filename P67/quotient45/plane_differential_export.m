// Prepared only: no native launch without root resource allocation.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
load "P67/data/quotient45/canonical_data.m";
load "P67/data/quotient45/jmap_q_data.m";
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
h:=Integers()!((GF(67)!PrimitiveRoot(67))^11);MS:=CuspidalSubspace(ModularSymbolsH(67,[h],2,1));decomp:=NewformDecomposition(MS);rows:=[];
for A in decomp do
 f:=qEigenform(A,7000);KF:=BaseRing(Parent(f));
 if KF cmpeq Rationals() then Append(~rows,[Rationals()!Coefficient(f,n):n in [0..6999]]);else
  dd:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(f,n)):n in [0..6999]];
  for k in [1..dd] do Append(~rows,[z[k]:z in cofs]);end for;
 end if;
end for;
Q0:=Matrix(Rationals(),rows);old:=ChangeRing(qcoeffs,Rationals());T:=Solution(Submatrix(Q0,1,1,45,320),old);full:=T*Q0;Q2:=ChangeRing(full,GF(2));assert Submatrix(Q2,1,1,45,1400) eq qcoeffs1400;
E,U:=EchelonForm(Submatrix(Q2,1,1,45,320));assert U*Submatrix(Q2,1,1,45,320) eq E;
LS<q>:=LaurentSeriesRing(GF(2),7000);F:=[LS!Eltseq(v)+O(q^7000):v in Rows(U*Q2)];ix:=[43,44,45];a,b,c:=Explode([F[i]:i in ix]);uu:=b/a;vv:=c/a;
if IsZero(Derivative(uu)) then assert not IsZero(Derivative(vv));ix:=[43,45,44];a,b,c:=Explode([F[i]:i in ix]);uu:=b/a;vv:=c/a;end if;
assert not IsZero(Derivative(uu));Dmax:=89-Valuation(a);assert Dmax le 46;
printf "NONIC67H3_PLANE_PROJECTION q_orders=%o max_plane_degree=%o separable=true time=%o\n",[Valuation(z):z in [a,b,c]],Dmax,Cputime(tt);
P<x,y,z>:=ProjectiveSpace(GF(2),2);R:=CoordinateRing(P);found:=false;
for d in [2..Dmax] do
 mons:=Setseq(MonomialsOfDegree(R,d));ser:=[Evaluate(m,[b,c,a]):m in mons];A:=Matrix(GF(2),[[Coefficient(s,n):n in [0..125*d]]:s in ser]);K:=Kernel(A);
 if Dimension(K) gt 0 then assert Dimension(K) eq 1;v:=K.1;peq:=&+[v[i]*mons[i]:i in [1..#mons]];assert IsIrreducible(peq);planeDeg:=d;found:=true;break;end if;
end for;assert found;
PrintFile("P67/data/quotient45/plane_projection_equation.m",Sprintf("N:=67;plane_degree:=%o;PlaneExponents:=%o;ProjectionRows:=Matrix(GF(2),3,45,%o);",planeDeg,[Exponents(m):m in Monomials(peq)],Eltseq(Matrix(GF(2),[Eltseq(U[i]):i in ix]))):Overwrite:=true);
AA<xx,yy>:=AffineSpace(GF(2),2);R2:=CoordinateRing(AA);CA:=Curve(AA,Evaluate(peq,[R2.1,R2.2,R2!1]));C:=ProjectiveClosure(CA);FF:=FunctionField(C);AF,am:=AlgorithmicFunctionField(FF);gx:=Genus(AF);assert gx eq 45;
uAF:=am(FF.1);vAF:=am(FF.2);KAF:=BaseRing(AF);printf "NONIC67H3_PLANE_AF_COORDINATES degree=%o ubase=%o vbase=%o ugen=%o vgen=%o time=%o\n",Degree(AF),uAF eq AF!KAF.1,vAF eq AF!KAF.1,uAF eq AF.1,vAF eq AF.1,Cputime(tt);if uAF eq AF!KAF.1 and vAF eq AF.1 then baseSeries:=uu;genSeries:=vv;baseCoord:=1;genCoord:=2;
elif vAF eq AF!KAF.1 and uAF eq AF.1 then baseSeries:=vv;genSeries:=uu;baseCoord:=2;genCoord:=1;
else error "Unexpected AF coordinate change";end if;
printf "NONIC67H3_PLANE_MODEL genus=%o plane_degree=%o AF_degree=%o time=%o\n",gx,planeDeg,Degree(AF),Cputime(tt);
Om,om:=SpaceOfHolomorphicDifferentials(AF);B:=[om(w):w in Basis(Om)];eta:=Differential(uAF);BC:=[w/eta:w in B];

function CodeK(w)
 return <[Integers()!a:a in Eltseq(Numerator(w))],[Integers()!a:a in Eltseq(Denominator(w))]>;
end function;
function CodeAF(w)return [CodeK(c):c in Eltseq(w)];end function;
out:="P67/data/quotient45/plane_differential_certified_data.m";
PrintFile(out,Sprintf("N:=67;genus:=45;plane_degree:=%o;AF_degree:=%o;AFBaseCoordinate:=%o;AFGeneratorCoordinate:=%o;PlaneExponents:=%o;ProjectionRows:=Matrix(GF(2),3,45,%o);ProjectionQSeries:=Matrix(GF(2),3,7000,%o);AFPolynomialData:=%o;DifferentialCoefficientData:=%o;",planeDeg,Degree(AF),baseCoord,genCoord,[Exponents(m):m in Monomials(peq)],Eltseq(Matrix(GF(2),[Eltseq(U[i]):i in ix])),Eltseq(Matrix(GF(2),[Eltseq((U*Q2)[i]):i in ix])),[CodeK(w):w in Eltseq(DefiningPolynomial(AF))],[CodeAF(w):w in BC]):Overwrite:=true);
printf "NONIC67H3_PLANE_DIFFERENTIAL_DATA_COMPLETE genus=45 plane_degree=%o AF_degree=%o basis_dimension=45 coefficient_base_coordinate=%o time=%o\n",planeDeg,Degree(AF),baseCoord,Cputime(tt);
PrintFile(out,"print \"NONIC67H3_PLANE_DIFFERENTIAL_DATA_READY\";");quit;
