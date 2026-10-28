// Prepared only: requires root allocation before any native launch.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
load "P67/data/quotient45/canonical_data.m";
N:=67;h:=Integers()!((GF(67)!PrimitiveRoot(67))^11);MS:=CuspidalSubspace(ModularSymbolsH(67,[h],2,1));decomp:=NewformDecomposition(MS);
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
Q0:=[];D0:=[];
for A in decomp do
 f:=qEigenform(A,1400);KF:=BaseRing(Parent(f));assert Coefficient(f,1) eq 1;chi:=(Coefficient(f,2)^2-Coefficient(f,4))/2;
 if KF cmpeq Rationals() then
  Append(~Q0,[Rationals()!Coefficient(f,n):n in [0..1399]]);Append(~D0,[Rationals()!(chi*Coefficient(f,n)):n in [0..1399]]);
 else
  d:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(f,n)):n in [0..1399]];dcofs:=[Flat(KF!(chi*Coefficient(f,n))):n in [0..1399]];
  for k in [1..d] do Append(~Q0,[z[k]:z in cofs]);Append(~D0,[z[k]:z in dcofs]);end for;
 end if;
end for;
Q0:=Matrix(Rationals(),Q0);D0:=Matrix(Rationals(),D0);assert Nrows(Q0) eq 45;
old:=ChangeRing(qcoeffs,Rationals());pre:=Submatrix(Q0,1,1,45,320);assert RowSpace(pre) eq RowSpace(old);T:=Solution(pre,old);assert T*pre eq old;
Qfull:=T*Q0;Dfull:=T*D0;assert &and[IsOdd(Denominator(a)):a in Eltseq(Qfull)];Q2:=ChangeRing(Qfull,GF(2));assert Submatrix(Q2,1,1,45,320) eq ChangeRing(qcoeffs,GF(2));
MdQ:=Solution(old,Submatrix(Dfull,1,1,45,320));assert MdQ*old eq Submatrix(Dfull,1,1,45,320);Md:=ChangeRing(MdQ,GF(2));I45:=IdentityMatrix(GF(2),45);assert Md^11 eq I45 and Md ne I45;
printf "NONIC67H3_JMAP_Q1400 same_basis=true diamond2_order=11 time=%o\n",Cputime(tt);
P:=ProjectiveSpace(GF(2),44);R:=CoordinateRing(P);mons2:=[R.(z[1])*R.(z[2]):z in pairs];eqs:=[&+[v[i]*mons2[i]:i in [1..1035]]:v in Rows(QuadRelations)];C:=Curve(P,eqs);IC:=Ideal(C);GB:=GroebnerBasis(IC);LM:=[LeadingMonomial(f):f in GB];
ls:=[&+[Md[i,j]*R.j:j in [1..45]]:i in [1..45]];assert &and[NormalForm(Evaluate(f,ls),IC) eq 0:f in eqs];
function StandardMonomials(d)
 ss:={R!1};
 for k in [1..d] do
  cand:={m*R.i:m in ss,i in [1..45]};ss:={m:m in cand|not exists{l:l in LM|IsDivisibleBy(m,l)}};
  assert #ss eq (k eq 1 select 45 else (2*k-1)*44);
 end for;
 return Sort(Setseq(ss));
end function;

S4:=StandardMonomials(4);S10:=StandardMonomials(10);assert #S4 eq 308 and #S10 eq 836;
PS<q>:=PowerSeriesRing(GF(2),1400);Q:=[PS!Eltseq(v)+O(q^1400):v in Rows(Q2)];
function ExpMat(ms)
 return Matrix(GF(2),[[Coefficient(f,n):n in [0..1399]] where f:=Evaluate(m,Q):m in ms]);
end function;
E4:=ExpMat(S4);E10:=ExpMat(S10);assert Rank(E4) eq 308 and Rank(E10) eq 836;
Delta:=q*&*[(1+q^(8*n))*(1+q^(16*n)):n in [1..174]]+O(q^1400);
DE4:=Matrix(GF(2),[[Coefficient(Delta*f,n):n in [0..1399]]where f:=PS!Eltseq(v)+O(q^1400):v in Rows(E4)]);
Z:=ZeroMatrix(GF(2),836,1400);KK:=Kernel(VerticalJoin(VerticalJoin(HorizontalJoin(E4,DE4),HorizontalJoin(E10,Z)),HorizontalJoin(Z,E10)));assert Dimension(KK) gt 0;
v:=KK.1;gv:=Vector(GF(2),Eltseq(v)[1..308]);av:=Vector(GF(2),Eltseq(v)[309..1144]);bv:=Vector(GF(2),Eltseq(v)[1145..1980]);assert gv ne 0 and gv*E4 eq av*E10 and gv*DE4 eq bv*E10;
num:=&+[av[i]*S10[i]:i in [1..836]];den:=&+[bv[i]*S10[i]:i in [1..836]];assert num ne 0 and den ne 0;
S2base:=StandardMonomials(2);assert #S2base eq 132;
out:="P67/data/quotient45/jmap_q_data.m";
PrintFile(out,Sprintf("N:=67;genus:=45;precision:=1400;SturmBoundWeight20:=1246;CanonicalJMapDegree:=10;AuxiliaryCanonicalDegree:=4;qcoeffs1400:=Matrix(GF(2),45,1400,%o);qcoeffs600:=Matrix(GF(2),45,600,%o);Diamond2Matrix:=Matrix(GF(2),45,45,%o);Std2Exponents:=%o;Std4Exponents:=%o;Std10Exponents:=%o;Gvec:=%o;NumeratorVec:=%o;DenominatorVec:=%o;",Eltseq(Q2),Eltseq(Submatrix(Q2,1,1,45,600)),Eltseq(Md),[Exponents(m):m in S2base],[Exponents(m):m in S4],[Exponents(m):m in S10],Eltseq(gv),Eltseq(av),Eltseq(bv)):Overwrite:=true);
printf "NONIC67H3_JMAP_IDENTITIES canonical_degree10 auxiliary_degree4 precision1400 Sturm1246 G_space_dimension=%o time=%o\n",Dimension(KK),Cputime(tt);
PrintFile(out,"print \"NONIC67H3_QUOTIENT_JMAP_Q_DATA_READY\";");print "NONIC67H3_JMAP_Q_EXPORT_COMPLETE";quit;
