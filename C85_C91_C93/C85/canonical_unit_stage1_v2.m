// Root-approved C85 geometry after exact completed finite-image character check.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "C85_C91_C93/data/character_probes/N85_winding_data.m";
assert Complete and 2 in FiniteDiamondUnits and 16 in FiniteDiamondUnits;
N:=85;g:=225;QPrec:=1300;TatePrec:=1800;
MS:=CuspidalSubspace(ModularSymbolsH(N,[1],2,1));assert Dimension(MS)eq g;decomp:=NewformDecomposition(MS);
printf "NONIC85_CANONICAL_UNIT_NEWFORMS dimensions=%o time=%o\n",[Dimension(A):A in decomp],Cputime(tt);
function Flat(z)
 if Type(z)eq FldRatElt or Type(z)eq RngIntElt then return[Rationals()|z];end if;return &cat[Flat(y):y in Eltseq(z)];
end function;
Q0:=[];D20:=[];D30:=[];
for A in decomp do
 lev:=Level(AssociatedNewSpace(A));assert lev in[17,85];first:=#Q0+1;ff:=qEigenform(A,QPrec);KF:=BaseRing(Parent(ff));chi2:=(Coefficient(ff,2)^2-Coefficient(ff,4))/2;chi3:=(Coefficient(ff,3)^2-Coefficient(ff,9))/3;
 if KF cmpeq Rationals()then
  Append(~Q0,[Rationals()!Coefficient(ff,n):n in[0..QPrec-1]]);Append(~D20,[Rationals()!(chi2*Coefficient(ff,n)):n in[0..QPrec-1]]);Append(~D30,[Rationals()!(chi3*Coefficient(ff,n)):n in[0..QPrec-1]]);
 else
  dd:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(ff,n)):n in[0..QPrec-1]];d2:=[Flat(KF!(chi2*Coefficient(ff,n))):n in[0..QPrec-1]];d3:=[Flat(KF!(chi3*Coefficient(ff,n))):n in[0..QPrec-1]];
  assert &and[#c eq dd:c in cofs];for k in[1..dd]do Append(~Q0,[c[k]:c in cofs]);Append(~D20,[c[k]:c in d2]);Append(~D30,[c[k]:c in d3]);end for;
 end if;
 if lev eq 17 then
  old:=Q0[first..#Q0];od2:=D20[first..#D20];od3:=D30[first..#D30];
  for rr in old do Append(~Q0,[n mod 5 eq 0 select rr[n div 5+1]else Rationals()!0:n in[0..QPrec-1]]);end for;
  for rr in od2 do Append(~D20,[n mod 5 eq 0 select rr[n div 5+1]else Rationals()!0:n in[0..QPrec-1]]);end for;
  for rr in od3 do Append(~D30,[n mod 5 eq 0 select rr[n div 5+1]else Rationals()!0:n in[0..QPrec-1]]);end for;
 end if;
 printf "NONIC85_CANONICAL_UNIT_EIGENFORM newlevel=%o dimension=%o accumulated_with_old_degeneracies=%o time=%o\n",lev,Dimension(A),#Q0,Cputime(tt);
end for;
QM0:=Matrix(Rationals(),Q0);DM20:=Matrix(Rationals(),D20);DM30:=Matrix(Rationals(),D30);assert Nrows(QM0)eq g;
den:=LCM([Denominator(a):a in Eltseq(QM0)]);QMZ:=Saturation(ChangeRing(den*QM0,Integers()));QMR:=ChangeRing(QMZ,Rationals());assert Nrows(QMZ)eq g and RowSpace(QMR)eq RowSpace(QM0);
Q2:=ChangeRing(QMZ,GF(2));assert Rank(Q2)eq g and Rank(Submatrix(Q2,1,1,g,577))eq g;TT:=Solution(QM0,QMR);assert TT*QM0 eq QMR;
Md2Q:=Solution(QMR,TT*DM20);Md3Q:=Solution(QMR,TT*DM30);assert Md2Q*QMR eq TT*DM20 and Md3Q*QMR eq TT*DM30;assert &and[IsOdd(Denominator(a)):a in Eltseq(Md2Q)cat Eltseq(Md3Q)];M2:=ChangeRing(Md2Q,GF(2));M3:=ChangeRing(Md3Q,GF(2));II:=IdentityMatrix(GF(2),g);assert M2^8 eq II and M2^4 ne II and M3^16 eq II and M3^8 ne II and M2*M3 eq M3*M2;
assert #{M2^a*M3^b:a in[0..1],b in[0..15]}eq 32;
printf "NONIC85_CANONICAL_UNIT_SATURATED_BASIS rank225=true precision1300=true Sturm_weight2=576 diamonds2_3_orders8_16=true grouporder32=true time=%o\n",Cputime(tt);
LS<q>:=LaurentSeriesRing(GF(2),TatePrec);xp:=q/(1+q)^2;yp:=q^2/(1+q)^3;s1:=LS!0;a6:=LS!0;
for n in[1..TatePrec div N+1]do
 if N*n+1 lt TatePrec then z:=q^(N*n+1);xp+:=z/(1+z)^2;yp+:=z^2/(1+z)^3;end if;
 if N*n-1 lt TatePrec then w:=q^(N*n-1);xp+:=w/(1+w)^2;yp+:=w/(1+w)^3;end if;
 if N*n lt TatePrec then z:=q^(N*n);s1+:=(GF(2)!n)*z/(1+z);a6+:=(GF(2)!((5*n^3+7*n^5)div 12))*z/(1+z);end if;
end for;
yp+:=s1;a4:=s1;assert Valuation(yp^2+xp*yp+xp^3+a4*xp+a6)ge TatePrec;
aa4:=xp^2+a4+yp;aaa2:=xp+(aa4/xp)^2+aa4/xp;b:=aaa2^3/xp^2;c:=(aaa2+xp)/xp;r:=b/c;s:=c^2/(b+c);t:=r*s+1;xu:=(s+r)/t;yu:=t/(s^2+s+r+1);assert AbsolutePrecision(xu)ge QPrec and AbsolutePrecision(yu)ge QPrec;
Q:=[LS!Eltseq(v)+O(q^QPrec):v in Rows(Q2)];low:=Min(0,Min(Valuation(xu),Valuation(yu))+1);RelPrec:=Min(QPrec,QPrec+Min(Valuation(xu),Valuation(yu)));assert RelPrec-1 gt 543;
function UnitRelations(uu)
 AA:=Matrix(GF(2),[[Coefficient(f,n):n in[low..RelPrec-1]]:f in Q cat[uu*f:f in Q]]);B:=BasisMatrix(Kernel(AA));assert Ncols(B)eq 2*g;
 for row in Rows(B)do assert Valuation(&+[row[i]*Q[i]+row[g+i]*uu*Q[i]:i in[1..g]])ge RelPrec;end for;return B;
end function;
RX:=UnitRelations(xu);RY:=UnitRelations(yu);assert Nrows(RX)ge 135 and Nrows(RY)ge 130;
// R2product rank is checked efficiently by the independent Python tensor step.
printf "NONIC85_CANONICAL_UNIT_RELATIONS x=%o y=%o exact_global_zero_order=%o section_degree_bounds=[538,543] product_R2_rank_pending_independent_tensor=true time=%o\n",Nrows(RX),Nrows(RY),RelPrec-1,Cputime(tt);
load "C85_C91_C93/data/N85_f9/records.m";
assert Complete and FieldDegree eq 9;PF:=PolynomialRing(GF(2));F<a>:=ext<GF(2)|PF!FieldPolynomial>;chosen:=false;
for rec in Records do
 E:=EllipticCurve([F!PF!v:v in rec[2]]);PP:=E![F!PF!v:v in rec[3]];assert Order(PP)eq 85;xy:=TateNormalForm_xy(E,PP);x0,y0:=Explode(xy);RR:=VerticalJoin(Matrix(F,[[F!v[i]+x0*v[g+i]:i in[1..g]]:v in Rows(RX)]),Matrix(F,[[F!v[i]+y0*v[g+i]:i in[1..g]]:v in Rows(RY)]));rk:=Rank(RR);printf "NONIC85_CANONICAL_UNIT_FINITE_RANK rank=%o expected224 time=%o\n",rk,Cputime(tt);if rk ne g-1 then continue;end if;
 zz:=Eltseq(Basis(Kernel(Transpose(RR)))[1]);ki:=Index([z ne 0:z in zz],true);zz:=[z/zz[ki]:z in zz];SampleE:=E;SampleP:=PP;SampleXY:=xy;SampleZ:=zz;chosen:=true;break;
end for;assert chosen;
function Code(z)return[Integers()!c:c in Eltseq(z,GF(2))];end function;
out:="C85_C91_C93/data/C85/canonical_unit_stage1_data.m";
PrintFile(out,Sprintf("N:=85;genus:=225;precision:=1300;UnitSeriesLow:=%o;UnitXSeries:=%o;UnitYSeries:=%o;qcoeffs:=Matrix(GF(2),225,1300,%o);Diamond2Matrix:=Matrix(GF(2),225,225,%o);Diamond3Matrix:=Matrix(GF(2),225,225,%o);UnitXRelations:=Matrix(GF(2),%o,450,%o);UnitYRelations:=Matrix(GF(2),%o,450,%o);SampleFieldPolynomial:=%o;SampleEllipticCoefficients:=%o;SampleLevelPoint:=%o;SampleUnitCoordinates:=%o;SampleCanonicalCoordinates:=%o;",low,[Integers()!Coefficient(xu,n):n in[low..QPrec-1]],[Integers()!Coefficient(yu,n):n in[low..QPrec-1]],Eltseq(Q2),Eltseq(M2),Eltseq(M3),Nrows(RX),Eltseq(RX),Nrows(RY),Eltseq(RY),[Integers()!c:c in Eltseq(DefiningPolynomial(F,GF(2)))],[Code(z):z in aInvariants(SampleE)],[Code(SampleP[i]):i in[1,2]],[Code(z):z in SampleXY],[Code(z):z in SampleZ]):Overwrite:=true);
PrintFile(out,"print \"NONIC85_CANONICAL_UNIT_STAGE1_DATA_READY\";");printf "NONIC85_CANONICAL_UNIT_STAGE1_COMPLETE generic_rank224=true no_function_field_normalization=true time=%o\n",Cputime(tt);quit;
