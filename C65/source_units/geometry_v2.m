// Root-approved C85 geometry after exact completed finite-image character check.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
assert assigned N and assigned OutDir;N:=StringToInteger(N);
g:=Dimension(CuspidalSubspace(ModularSymbolsH(N,[1],2,1)));
QPrec:=N eq 65 select 1100 else (N eq 117 select 2400 else 1800);TatePrec:=QPrec+500;
MS:=CuspidalSubspace(ModularSymbolsH(N,[1],2,1));decomp:=NewformDecomposition(MS);
D1:=2;D2:=N mod 3 eq 0 select 5 else 3;
printf "SOURCE_UNIT_START level=%o genus=%o precision=%o diamonds=%o time=%o\n",N,g,QPrec,[D1,D2],Cputime(tt);
function Flat(z)
 if Type(z)eq FldRatElt or Type(z)eq RngIntElt then return[Rationals()|z];end if;return &cat[Flat(y):y in Eltseq(z)];
end function;
Q0:=[];D20:=[];D30:=[];
for A in decomp do
 lev:=Level(AssociatedNewSpace(A));assert N mod lev eq 0;first:=#Q0+1;ff:=qEigenform(A,QPrec);KF:=BaseRing(Parent(ff));chi2:=(Coefficient(ff,D1)^2-Coefficient(ff,D1^2))/D1;chi3:=(Coefficient(ff,D2)^2-Coefficient(ff,D2^2))/D2;
 if KF cmpeq Rationals()then
  Append(~Q0,[Rationals()!Coefficient(ff,n):n in[0..QPrec-1]]);Append(~D20,[Rationals()!(chi2*Coefficient(ff,n)):n in[0..QPrec-1]]);Append(~D30,[Rationals()!(chi3*Coefficient(ff,n)):n in[0..QPrec-1]]);
 else
  dd:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(ff,n)):n in[0..QPrec-1]];d2:=[Flat(KF!(chi2*Coefficient(ff,n))):n in[0..QPrec-1]];d3:=[Flat(KF!(chi3*Coefficient(ff,n))):n in[0..QPrec-1]];
  assert &and[#c eq dd:c in cofs];for k in[1..dd]do Append(~Q0,[c[k]:c in cofs]);Append(~D20,[c[k]:c in d2]);Append(~D30,[c[k]:c in d3]);end for;
 end if;
 if lev ne N then
  old:=Q0[first..#Q0];od2:=D20[first..#D20];od3:=D30[first..#D30];
  for dv in Divisors(N div lev) do
   if dv eq 1 then continue;end if;
   for rr in old do Append(~Q0,[n mod dv eq 0 select rr[n div dv+1]else Rationals()!0:n in[0..QPrec-1]]);end for;
   for rr in od2 do Append(~D20,[n mod dv eq 0 select rr[n div dv+1]else Rationals()!0:n in[0..QPrec-1]]);end for;
   for rr in od3 do Append(~D30,[n mod dv eq 0 select rr[n div dv+1]else Rationals()!0:n in[0..QPrec-1]]);end for;
  end for;
 end if;
 printf "NONIC85_CANONICAL_UNIT_EIGENFORM newlevel=%o dimension=%o accumulated_with_old_degeneracies=%o time=%o\n",lev,Dimension(A),#Q0,Cputime(tt);
end for;
QM0:=Matrix(Rationals(),Q0);DM20:=Matrix(Rationals(),D20);DM30:=Matrix(Rationals(),D30);assert Nrows(QM0)eq g;
den:=LCM([Denominator(a):a in Eltseq(QM0)]);QMZ:=Saturation(ChangeRing(den*QM0,Integers()));QMR:=ChangeRing(QMZ,Rationals());assert Nrows(QMZ)eq g and RowSpace(QMR)eq RowSpace(QM0);
Q2:=ChangeRing(QMZ,GF(2));assert Rank(Q2)eq g;TT:=Solution(QM0,QMR);assert TT*QM0 eq QMR;
Md2Q:=Solution(QMR,TT*DM20);Md3Q:=Solution(QMR,TT*DM30);assert Md2Q*QMR eq TT*DM20 and Md3Q*QMR eq TT*DM30;assert &and[IsOdd(Denominator(a)):a in Eltseq(Md2Q)cat Eltseq(Md3Q)];M2:=ChangeRing(Md2Q,GF(2));M3:=ChangeRing(Md3Q,GF(2));II:=IdentityMatrix(GF(2),g);assert M2*M3 eq M3*M2;
printf "SOURCE_UNIT_SATURATION level=%o genus=%o rank=%o exact_Q_span=true time=%o\n",N,g,Rank(Q2),Cputime(tt);
LS<q>:=LaurentSeriesRing(GF(2),TatePrec);xp:=q/(1+q)^2;yp:=q^2/(1+q)^3;s1:=LS!0;a6:=LS!0;
for n in[1..TatePrec div N+1]do
 if N*n+1 lt TatePrec then z:=q^(N*n+1);xp+:=z/(1+z)^2;yp+:=z^2/(1+z)^3;end if;
 if N*n-1 lt TatePrec then w:=q^(N*n-1);xp+:=w/(1+w)^2;yp+:=w/(1+w)^3;end if;
 if N*n lt TatePrec then z:=q^(N*n);s1+:=(GF(2)!n)*z/(1+z);a6+:=(GF(2)!((5*n^3+7*n^5)div 12))*z/(1+z);end if;
end for;
yp+:=s1;a4:=s1;assert Valuation(yp^2+xp*yp+xp^3+a4*xp+a6)ge TatePrec;
aa4:=xp^2+a4+yp;aaa2:=xp+(aa4/xp)^2+aa4/xp;b:=aaa2^3/xp^2;c:=(aaa2+xp)/xp;r:=b/c;s:=c^2/(b+c);t:=r*s+1;xu:=(s+r)/t;yu:=t/(s^2+s+r+1);assert AbsolutePrecision(xu)ge QPrec and AbsolutePrecision(yu)ge QPrec;
Q:=[LS!Eltseq(v)+O(q^QPrec):v in Rows(Q2)];low:=Min(0,Min(Valuation(xu),Valuation(yu))+1);RelPrec:=Min(QPrec,QPrec+Min(Valuation(xu),Valuation(yu)));assert RelPrec gt 0;
function UnitRelations(uu)
 AA:=Matrix(GF(2),[[Coefficient(f,n):n in[low..RelPrec-1]]:f in Q cat[uu*f:f in Q]]);B:=BasisMatrix(Kernel(AA));assert Ncols(B)eq 2*g;
 for row in Rows(B)do assert Valuation(&+[row[i]*Q[i]+row[g+i]*uu*Q[i]:i in[1..g]])ge RelPrec;end for;return B;
end function;
RX:=UnitRelations(xu);RY:=UnitRelations(yu);
printf "SOURCE_UNIT_RELATIONS level=%o genus=%o x=%o y=%o vanishing_precision=%o time=%o\n",N,g,Nrows(RX),Nrows(RY),RelPrec,Cputime(tt);
function Code(z)return[Integers()!c:c in Eltseq(z,GF(2))];end function;
load "C65/data/census/N65_f9_records.m";
PF:=PolynomialRing(GF(2));F<a>:=ext<GF(2)|PF!FieldPolynomial>;chosen:=false;
for rec in Records do
 E:=EllipticCurve([F!PF!v:v in rec[2]]);PP:=E![F!PF!v:v in rec[3]];assert Order(PP)eq N;
 xy:=TateNormalForm_xy(E,PP);x0,y0:=Explode(xy);
 RR:=VerticalJoin(Matrix(F,[[F!v[i]+x0*v[g+i]:i in[1..g]]:v in Rows(RX)]),Matrix(F,[[F!v[i]+y0*v[g+i]:i in[1..g]]:v in Rows(RY)]));rk:=Rank(RR);
 printf "SOURCE_UNIT_SAMPLE level=%o rank=%o expected=%o time=%o\n",N,rk,g-1,Cputime(tt);
 if rk eq g-1 then chosen:=true;break;end if;
end for;
out:=OutDir cat "/canonical_unit_data.m";
PrintFile(out,Sprintf("N:=%o;genus:=%o;precision:=%o;UnitSeriesLow:=%o;UnitXSeries:=%o;UnitYSeries:=%o;qcoeffs:=Matrix(GF(2),%o,%o,%o);DiamondUnits:=%o;Diamond1Matrix:=Matrix(GF(2),%o,%o,%o);Diamond2Matrix:=Matrix(GF(2),%o,%o,%o);UnitXRelations:=Matrix(GF(2),%o,%o,%o);UnitYRelations:=Matrix(GF(2),%o,%o,%o);GenericRecoverySamplePassed:=%o;",N,g,QPrec,low,[Integers()!Coefficient(xu,n):n in[low..QPrec-1]],[Integers()!Coefficient(yu,n):n in[low..QPrec-1]],g,QPrec,Eltseq(Q2),[D1,D2],g,g,Eltseq(M2),g,g,Eltseq(M3),Nrows(RX),2*g,Eltseq(RX),Nrows(RY),2*g,Eltseq(RY),chosen):Overwrite:=true);
printf "SOURCE_UNIT_COMPLETE level=%o genus=%o sample_recovery=%o time=%o\n",N,g,chosen,Cputime(tt);quit;
