// Root-allocated first canonical/unit stage only; no function-field normalization.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
N:=83;g:=247;QPrec:=1300;TatePrec:=1800;
MS:=CuspidalSubspace(ModularSymbolsH(N,[1],2,1));assert Dimension(MS)eq g;decomp:=NewformDecomposition(MS);
printf "NONIC83_CANONICAL_UNIT_NEWFORMS dimensions=%o time=%o\n",[Dimension(A):A in decomp],Cputime(tt);
function Flat(z)
 if Type(z)eq FldRatElt or Type(z)eq RngIntElt then return [Rationals()|z];end if;return &cat[Flat(y):y in Eltseq(z)];
end function;
Q0:=[];D0:=[];
for A in decomp do
 ff:=qEigenform(A,QPrec);KF:=BaseRing(Parent(ff));chi:=(Coefficient(ff,2)^2-Coefficient(ff,4))/2;
 if KF cmpeq Rationals()then
  Append(~Q0,[Rationals()!Coefficient(ff,n):n in[0..QPrec-1]]);Append(~D0,[Rationals()!(chi*Coefficient(ff,n)):n in[0..QPrec-1]]);
 else
  dd:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(ff,n)):n in[0..QPrec-1]];dcofs:=[Flat(KF!(chi*Coefficient(ff,n))):n in[0..QPrec-1]];
  assert &and[#c eq dd:c in cofs]and &and[#c eq dd:c in dcofs];for k in[1..dd]do Append(~Q0,[c[k]:c in cofs]);Append(~D0,[c[k]:c in dcofs]);end for;
 end if;
 printf "NONIC83_CANONICAL_UNIT_EIGENFORM dimension=%o accumulated=%o time=%o\n",Dimension(A),#Q0,Cputime(tt);
end for;
QM0:=Matrix(Rationals(),Q0);DM0:=Matrix(Rationals(),D0);assert Nrows(QM0)eq g;
den:=LCM([Denominator(a):a in Eltseq(QM0)]);QMZ:=Saturation(ChangeRing(den*QM0,Integers()));QMR:=ChangeRing(QMZ,Rationals());assert Nrows(QMZ)eq g and RowSpace(QMR)eq RowSpace(QM0);
Q2:=ChangeRing(QMZ,GF(2));assert Rank(Q2)eq g and Rank(Submatrix(Q2,1,1,g,575))eq g;
TT:=Solution(QM0,QMR);assert TT*QM0 eq QMR;Dfull:=TT*DM0;MdQ:=Solution(QMR,Dfull);assert MdQ*QMR eq Dfull and &and[IsOdd(Denominator(a)):a in Eltseq(MdQ)];Md:=ChangeRing(MdQ,GF(2));assert Md^41 eq IdentityMatrix(GF(2),g)and Md ne IdentityMatrix(GF(2),g);
printf "NONIC83_CANONICAL_UNIT_SATURATED_BASIS rank247=true precision1300=true Sturm_weight2=574 diamond2order41=true time=%o\n",Cputime(tt);
LS<q>:=LaurentSeriesRing(GF(2),TatePrec);xp:=q/(1+q)^2;yp:=q^2/(1+q)^3;s1:=LS!0;a6:=LS!0;
for n in[1..TatePrec div N+1]do
 if N*n+1 lt TatePrec then z:=q^(N*n+1);xp+:=z/(1+z)^2;yp+:=z^2/(1+z)^3;end if;
 if N*n-1 lt TatePrec then w:=q^(N*n-1);xp+:=w/(1+w)^2;yp+:=w/(1+w)^3;end if;
 if N*n lt TatePrec then z:=q^(N*n);s1+:=(GF(2)!n)*z/(1+z);a6+:=(GF(2)!((5*n^3+7*n^5)div 12))*z/(1+z);end if;
end for;
yp+:=s1;a4:=s1;assert Valuation(yp^2+xp*yp+xp^3+a4*xp+a6)ge TatePrec;
// The universal Tate point u=q on E(q^83) has exact order83.
aa4:=xp^2+a4+yp;aaa2:=xp+(aa4/xp)^2+aa4/xp;b:=aaa2^3/xp^2;c:=(aaa2+xp)/xp;r:=b/c;s:=c^2/(b+c);t:=r*s+1;xu:=(s+r)/t;yu:=t/(s^2+s+r+1);
assert AbsolutePrecision(xu)ge QPrec and AbsolutePrecision(yu)ge QPrec and Valuation(xu)ge -95 and Valuation(yu)ge -95;
printf "NONIC83_TATE_UNIT_SERIES valuation_x=%o valuation_y=%o absolute_precision_x=%o absolute_precision_y=%o time=%o\n",Valuation(xu),Valuation(yu),AbsolutePrecision(xu),AbsolutePrecision(yu),Cputime(tt);
Q:=[LS!Eltseq(v)+O(q^QPrec):v in Rows(Q2)];low:=Min(0,Min(Valuation(xu),Valuation(yu))+1);RelPrec:=Min(QPrec,QPrec+Min(Valuation(xu),Valuation(yu)));assert RelPrec-1 gt 587;
function UnitRelations(uu)
 AA:=Matrix(GF(2),[[Coefficient(f,n):n in[low..RelPrec-1]]:f in Q cat[uu*f:f in Q]]);B:=BasisMatrix(Kernel(AA));assert Ncols(B)eq 2*g;
 for row in Rows(B)do assert Valuation(&+[row[i]*Q[i]+row[g+i]*uu*Q[i]:i in[1..g]])ge RelPrec;end for;
 return B;
end function;
RX:=UnitRelations(xu);RY:=UnitRelations(yu);assert Nrows(RX)ge 157 and Nrows(RY)ge 152;
printf "NONIC83_CANONICAL_UNIT_RELATIONS x=%o y=%o exact_global_order_lower_bound=%o section_degree_bounds=[582,587] time=%o\n",Nrows(RX),Nrows(RY),RelPrec-1,Cputime(tt);
F<a>:=GF(2^9);chosen:=false;
for jj in F do
 for E in Twists(EllipticCurveFromjInvariant(jj))do
  if #E mod 83 ne 0 then continue;end if;GG,mp:=AbelianGroup(E);iv:=Invariants(GG);m:=iv[#iv];if m mod 83 ne 0 then continue;end if;PP:=(m div 83)*mp(GG.#iv);assert Order(PP)eq 83;
  xy:=TateNormalForm_xy(E,PP);x0,y0:=Explode(xy);RR:=VerticalJoin(Matrix(F,[[F!v[i]+x0*v[g+i]:i in[1..g]]:v in Rows(RX)]),Matrix(F,[[F!v[i]+y0*v[g+i]:i in[1..g]]:v in Rows(RY)]));rk:=Rank(RR);
  printf "NONIC83_CANONICAL_UNIT_FINITE_RANK rank=%o expected246 time=%o\n",rk,Cputime(tt);if rk ne g-1 then continue;end if;
  zz:=Eltseq(Basis(Kernel(Transpose(RR)))[1]);ki:=Index([z ne 0:z in zz],true);zz:=[z/zz[ki]:z in zz];SampleE:=E;SampleP:=PP;SampleXY:=xy;SampleZ:=zz;chosen:=true;break;
 end for;if chosen then break;end if;
end for;assert chosen;
out:="P83/data/canonical_units/canonical_unit_stage1_data.m";
function Code(z)return [Integers()!c:c in Eltseq(z,GF(2))];end function;
PrintFile(out,Sprintf("N:=83;genus:=247;precision:=1300;UnitSeriesLow:=%o;UnitXSeries:=%o;UnitYSeries:=%o;qcoeffs:=Matrix(GF(2),247,1300,%o);Diamond2Matrix:=Matrix(GF(2),247,247,%o);UnitXRelations:=Matrix(GF(2),%o,494,%o);UnitYRelations:=Matrix(GF(2),%o,494,%o);SampleFieldPolynomial:=%o;SampleEllipticCoefficients:=%o;SampleLevelPoint:=%o;SampleUnitCoordinates:=%o;SampleCanonicalCoordinates:=%o;",low,[Integers()!Coefficient(xu,n):n in[low..QPrec-1]],[Integers()!Coefficient(yu,n):n in[low..QPrec-1]],Eltseq(Q2),Eltseq(Md),Nrows(RX),Eltseq(RX),Nrows(RY),Eltseq(RY),[Integers()!c:c in Eltseq(DefiningPolynomial(F,GF(2)))],[Code(z):z in aInvariants(SampleE)],[Code(SampleP[i]):i in[1,2]],[Code(z):z in SampleXY],[Code(z):z in SampleZ]):Overwrite:=true);
PrintFile(out,"print \"NONIC83_CANONICAL_UNIT_STAGE1_DATA_READY\";");printf "NONIC83_CANONICAL_UNIT_STAGE1_COMPLETE generic_rank246=true no_function_field_normalization=true time=%o\n",Cputime(tt);quit;
