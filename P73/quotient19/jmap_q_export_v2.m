SetColumns(0);SetQuitOnError(true);SetMemoryLimit(48*10^9);tt:=Cputime();
load "P73/data/quotient19/canonical_data.m";
N:=73;h:=Integers()!((GF(73)!PrimitiveRoot(73))^4);MS:=CuspidalSubspace(ModularSymbolsH(73,[h],2,1));decomp:=NewformDecomposition(MS);
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
Q0:=[];D0:=[];
for A in decomp do
 f:=qEigenform(A,600);KF:=BaseRing(Parent(f));assert Coefficient(f,1) eq 1;chi:=(Coefficient(f,5)^2-Coefficient(f,25))/5;
 if KF cmpeq Rationals() then
  Append(~Q0,[Rationals()!Coefficient(f,n):n in [0..599]]);Append(~D0,[Rationals()!(chi*Coefficient(f,n)):n in [0..599]]);
 else
  d:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(f,n)):n in [0..599]];dcofs:=[Flat(KF!(chi*Coefficient(f,n))):n in [0..599]];
  for k in [1..d] do Append(~Q0,[z[k]:z in cofs]);Append(~D0,[z[k]:z in dcofs]);end for;
 end if;
end for;
Q0:=Matrix(Rationals(),Q0);D0:=Matrix(Rationals(),D0);assert Nrows(Q0) eq 19;
old:=ChangeRing(qcoeffs,Rationals());pre:=Submatrix(Q0,1,1,19,130);assert RowSpace(pre) eq RowSpace(old);T:=Solution(pre,old);assert T*pre eq old;
Qfull:=T*Q0;Dfull:=T*D0;assert &and[IsOdd(Denominator(a)):a in Eltseq(Qfull)];Q2:=ChangeRing(Qfull,GF(2));assert Submatrix(Q2,1,1,19,130) eq ChangeRing(qcoeffs,GF(2));
MdQ:=Solution(old,Submatrix(Dfull,1,1,19,130));assert MdQ*old eq Submatrix(Dfull,1,1,19,130);Md:=ChangeRing(MdQ,GF(2));I27:=IdentityMatrix(GF(2),19);assert Md^4 eq I27 and Md^2 ne I27;
printf "NONIC73_JMAP_Q600 same_basis=true diamond5_order=4 time=%o\n",Cputime(tt);
P:=ProjectiveSpace(GF(2),18);R:=CoordinateRing(P);mons2:=[R.(z[1])*R.(z[2]):z in pairs];eqs:=[&+[v[i]*mons2[i]:i in [1..190]]:v in Rows(QuadRelations)];C:=Curve(P,eqs);IC:=Ideal(C);GB:=GroebnerBasis(IC);LM:=[LeadingMonomial(f):f in GB];
ls:=[&+[Md[i,j]*R.j:j in [1..19]]:i in [1..19]];assert &and[NormalForm(Evaluate(f,ls),IC) eq 0:f in eqs];
function StandardMonomials(d)
 ss:={R!1};
 for k in [1..d] do
  cand:={m*R.i:m in ss,i in [1..19]};ss:={m:m in cand|not exists{l:l in LM|IsDivisibleBy(m,l)}};
  assert #ss eq (k eq 1 select 19 else (2*k-1)*18);
 end for;
 return Sort(Setseq(ss));
end function;

S3:=StandardMonomials(3);S9:=StandardMonomials(9);assert #S3 eq 90 and #S9 eq 306;
PS<q>:=PowerSeriesRing(GF(2),600);Q:=[PS!Eltseq(v)+O(q^600):v in Rows(Q2)];
function ExpMat(ms)
 return Matrix(GF(2),[[Coefficient(f,n):n in [0..599]] where f:=Evaluate(m,Q):m in ms]);
end function;
E3:=ExpMat(S3);E9:=ExpMat(S9);assert Rank(E3) eq 90 and Rank(E9) eq 306;
Delta:=q*&*[(1+q^(8*n))*(1+q^(16*n)):n in [1..74]]+O(q^600);
DE3:=Matrix(GF(2),[[Coefficient(Delta*f,n):n in [0..599]]where f:=PS!Eltseq(v)+O(q^600):v in Rows(E3)]);
Z:=ZeroMatrix(GF(2),306,600);KK:=Kernel(VerticalJoin(VerticalJoin(HorizontalJoin(E3,DE3),HorizontalJoin(E9,Z)),HorizontalJoin(Z,E9)));assert Dimension(KK) gt 0;
v:=KK.1;gv:=Vector(GF(2),Eltseq(v)[1..90]);av:=Vector(GF(2),Eltseq(v)[91..396]);bv:=Vector(GF(2),Eltseq(v)[397..702]);assert gv ne 0 and gv*E3 eq av*E9 and gv*DE3 eq bv*E9;
num:=&+[av[i]*S9[i]:i in [1..306]];den:=&+[bv[i]*S9[i]:i in [1..306]];assert num ne 0 and den ne 0;
out:="P73/data/quotient19/jmap_q_data.m";
PrintFile(out,Sprintf("N:=73;genus:=19;precision:=600;SturmBoundWeight18:=444;CanonicalJMapDegree:=9;AuxiliaryCanonicalDegree:=3;qcoeffs600:=Matrix(GF(2),19,600,%o);Diamond5Matrix:=Matrix(GF(2),19,19,%o);Std3Exponents:=%o;Std9Exponents:=%o;Gvec:=%o;NumeratorVec:=%o;DenominatorVec:=%o;",Eltseq(Q2),Eltseq(Md),[Exponents(m):m in S3],[Exponents(m):m in S9],Eltseq(gv),Eltseq(av),Eltseq(bv)):Overwrite:=true);
printf "NONIC73_JMAP_IDENTITIES canonical_degree=9 numerator=E4cubed_times_G denominator=Delta_times_G simultaneous=true precision=600 Sturm_bound=444 G_space_dimension=%o time=%o\n",Dimension(KK),Cputime(tt);
PrintFile(out,"print \"NONIC73_QUOTIENT_JMAP_Q_DATA_READY\";");print "NONIC73_JMAP_Q_EXPORT_COMPLETE";quit;
