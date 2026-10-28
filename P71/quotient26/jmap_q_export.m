// Prepared only: requires root allocation before any native launch.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
load "P71/data/quotient26/canonical_data.m";
N:=71;h:=Integers()!((GF(71)!PrimitiveRoot(71))^5);MS:=CuspidalSubspace(ModularSymbolsH(71,[h],2,1));decomp:=NewformDecomposition(MS);
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
Q0:=[];D0:=[];
for A in decomp do
 f:=qEigenform(A,600);KF:=BaseRing(Parent(f));assert Coefficient(f,1) eq 1;chi:=(Coefficient(f,7)^2-Coefficient(f,49))/7;
 if KF cmpeq Rationals() then
  Append(~Q0,[Rationals()!Coefficient(f,n):n in [0..599]]);Append(~D0,[Rationals()!(chi*Coefficient(f,n)):n in [0..599]]);
 else
  d:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(f,n)):n in [0..599]];dcofs:=[Flat(KF!(chi*Coefficient(f,n))):n in [0..599]];
  for k in [1..d] do Append(~Q0,[z[k]:z in cofs]);Append(~D0,[z[k]:z in dcofs]);end for;
 end if;
end for;
Q0:=Matrix(Rationals(),Q0);D0:=Matrix(Rationals(),D0);assert Nrows(Q0) eq 26;
old:=ChangeRing(qcoeffs,Rationals());pre:=Submatrix(Q0,1,1,26,130);assert RowSpace(pre) eq RowSpace(old);T:=Solution(pre,old);assert T*pre eq old;
Qfull:=T*Q0;Dfull:=T*D0;assert &and[IsOdd(Denominator(a)):a in Eltseq(Qfull)];Q2:=ChangeRing(Qfull,GF(2));assert Submatrix(Q2,1,1,26,130) eq ChangeRing(qcoeffs,GF(2));
MdQ:=Solution(old,Submatrix(Dfull,1,1,26,130));assert MdQ*old eq Submatrix(Dfull,1,1,26,130);Md:=ChangeRing(MdQ,GF(2));I26:=IdentityMatrix(GF(2),26);assert Md^5 eq I26 and Md ne I26;
printf "NONIC71_JMAP_Q600 same_basis=true diamond7_order=5 time=%o\n",Cputime(tt);
P:=ProjectiveSpace(GF(2),25);R:=CoordinateRing(P);mons2:=[R.(z[1])*R.(z[2]):z in pairs];eqs:=[&+[v[i]*mons2[i]:i in [1..351]]:v in Rows(QuadRelations)];C:=Curve(P,eqs);IC:=Ideal(C);GB:=GroebnerBasis(IC);LM:=[LeadingMonomial(f):f in GB];
ls:=[&+[Md[i,j]*R.j:j in [1..26]]:i in [1..26]];assert &and[NormalForm(Evaluate(f,ls),IC) eq 0:f in eqs];
function StandardMonomials(d)
 ss:={R!1};
 for k in [1..d] do
  cand:={m*R.i:m in ss,i in [1..26]};ss:={m:m in cand|not exists{l:l in LM|IsDivisibleBy(m,l)}};
  assert #ss eq (k eq 1 select 26 else (2*k-1)*25);
 end for;
 return Sort(Setseq(ss));
end function;

S2:=StandardMonomials(2);S8:=StandardMonomials(8);assert #S2 eq 75 and #S8 eq 375;
PS<q>:=PowerSeriesRing(GF(2),600);Q:=[PS!Eltseq(v)+O(q^600):v in Rows(Q2)];
function ExpMat(ms)
 return Matrix(GF(2),[[Coefficient(f,n):n in [0..599]] where f:=Evaluate(m,Q):m in ms]);
end function;
E2:=ExpMat(S2);E8:=ExpMat(S8);assert Rank(E2) eq 75 and Rank(E8) eq 375;
Delta:=q*&*[(1+q^(8*n))*(1+q^(16*n)):n in [1..74]]+O(q^600);
DE2:=Matrix(GF(2),[[Coefficient(Delta*f,n):n in [0..599]]where f:=PS!Eltseq(v)+O(q^600):v in Rows(E2)]);
Z:=ZeroMatrix(GF(2),375,600);KK:=Kernel(VerticalJoin(VerticalJoin(HorizontalJoin(E2,DE2),HorizontalJoin(E8,Z)),HorizontalJoin(Z,E8)));assert Dimension(KK) gt 0;
v:=KK.1;gv:=Vector(GF(2),Eltseq(v)[1..75]);av:=Vector(GF(2),Eltseq(v)[76..450]);bv:=Vector(GF(2),Eltseq(v)[451..825]);assert gv ne 0 and gv*E2 eq av*E8 and gv*DE2 eq bv*E8;
num:=&+[av[i]*S8[i]:i in [1..375]];den:=&+[bv[i]*S8[i]:i in [1..375]];assert num ne 0 and den ne 0;
out:="P71/data/quotient26/jmap_q_data.m";
PrintFile(out,Sprintf("N:=71;genus:=26;precision:=600;SturmBoundWeight16:=480;CanonicalJMapDegree:=8;AuxiliaryCanonicalDegree:=2;qcoeffs600:=Matrix(GF(2),26,600,%o);Diamond7Matrix:=Matrix(GF(2),26,26,%o);Std2Exponents:=%o;Std8Exponents:=%o;Gvec:=%o;NumeratorVec:=%o;DenominatorVec:=%o;",Eltseq(Q2),Eltseq(Md),[Exponents(m):m in S2],[Exponents(m):m in S8],Eltseq(gv),Eltseq(av),Eltseq(bv)):Overwrite:=true);
printf "NONIC71_JMAP_IDENTITIES canonical_degree=8 numerator=E4cubed_times_G denominator=Delta_times_G simultaneous=true precision=600 Sturm_bound=480 G_space_dimension=%o time=%o\n",Dimension(KK),Cputime(tt);
PrintFile(out,"print \"NONIC71_QUOTIENT_JMAP_Q_DATA_READY\";");print "NONIC71_JMAP_Q_EXPORT_COMPLETE";quit;
