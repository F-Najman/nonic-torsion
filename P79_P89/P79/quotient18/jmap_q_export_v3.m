SetColumns(0);SetQuitOnError(true);SetMemoryLimit(48*10^9);tt:=Cputime();
load "P79_P89/data/P79/quotient18/canonical_data.m";
N:=79;h:=Integers()!((GF(79)!PrimitiveRoot(79))^3);MS:=CuspidalSubspace(ModularSymbolsH(79,[h],2,1));decomp:=NewformDecomposition(MS);
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
Q0:=[];D0:=[];
for A in decomp do
 f:=qEigenform(A,600);KF:=BaseRing(Parent(f));assert Coefficient(f,1) eq 1;chi:=(Coefficient(f,3)^2-Coefficient(f,9))/3;
 if KF cmpeq Rationals() then
  Append(~Q0,[Rationals()!Coefficient(f,n):n in [0..599]]);Append(~D0,[Rationals()!(chi*Coefficient(f,n)):n in [0..599]]);
 else
  d:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(f,n)):n in [0..599]];dcofs:=[Flat(KF!(chi*Coefficient(f,n))):n in [0..599]];
  for k in [1..d] do Append(~Q0,[z[k]:z in cofs]);Append(~D0,[z[k]:z in dcofs]);end for;
 end if;
end for;
Q0:=Matrix(Rationals(),Q0);D0:=Matrix(Rationals(),D0);assert Nrows(Q0) eq 18;
old:=ChangeRing(qcoeffs,Rationals());pre:=Submatrix(Q0,1,1,18,130);assert RowSpace(pre) eq RowSpace(old);T:=Solution(pre,old);assert T*pre eq old;
Qfull:=T*Q0;Dfull:=T*D0;assert &and[IsOdd(Denominator(a)):a in Eltseq(Qfull)];Q2:=ChangeRing(Qfull,GF(2));assert Submatrix(Q2,1,1,18,130) eq ChangeRing(qcoeffs,GF(2));
MdQ:=Solution(old,Submatrix(Dfull,1,1,18,130));assert MdQ*old eq Submatrix(Dfull,1,1,18,130);Md:=ChangeRing(MdQ,GF(2));I27:=IdentityMatrix(GF(2),18);assert Md^3 eq I27 and Md ne I27;
printf "NONIC79_JMAP_Q600 same_basis=true diamond3_order=3 time=%o\n",Cputime(tt);
P:=ProjectiveSpace(GF(2),17);R:=CoordinateRing(P);mons2:=[R.(z[1])*R.(z[2]):z in pairs];eqs:=[&+[v[i]*mons2[i]:i in [1..171]]:v in Rows(QuadRelations)];C:=Curve(P,eqs);IC:=Ideal(C);GB:=GroebnerBasis(IC);LM:=[LeadingMonomial(f):f in GB];
ls:=[&+[Md[i,j]*R.j:j in [1..18]]:i in [1..18]];assert &and[NormalForm(Evaluate(f,ls),IC) eq 0:f in eqs];
function StandardMonomials(d)
 ss:={R!1};
 for k in [1..d] do
  cand:={m*R.i:m in ss,i in [1..18]};ss:={m:m in cand|not exists{l:l in LM|IsDivisibleBy(m,l)}};
  assert #ss eq (k eq 1 select 18 else (2*k-1)*17);
 end for;
 return Sort(Setseq(ss));
end function;
S2:=StandardMonomials(2);S8:=StandardMonomials(8);assert #S2 eq 51 and #S8 eq 255;
PS<q>:=PowerSeriesRing(GF(2),600);Q:=[PS!Eltseq(v)+O(q^600):v in Rows(Q2)];
function ExpMat(ms)
 return Matrix(GF(2),[[Coefficient(f,n):n in [0..599]] where f:=Evaluate(m,Q):m in ms]);
end function;
E2:=ExpMat(S2);E8:=ExpMat(S8);assert Rank(E2) eq 51 and Rank(E8) eq 255;
K:=Kernel(VerticalJoin(E2,E8));assert Dimension(K) gt 0;
v:=K.1;gv:=Vector(GF(2),Eltseq(v)[1..51]);av:=Vector(GF(2),Eltseq(v)[52..306]);assert gv ne 0 and gv*E2 eq av*E8;
G:=PS!Eltseq(gv*E2)+O(q^600);Delta:=q*&*[(1+q^(8*n))*(1+q^(16*n)):n in [1..74]]+O(q^600);
bv:=Solution(E8,Vector(GF(2),[Coefficient(Delta*G,n):n in [0..599]]));assert bv*E8 eq Vector(GF(2),[Coefficient(Delta*G,n):n in [0..599]]);
num:=&+[av[i]*S8[i]:i in [1..255]];den:=&+[bv[i]*S8[i]:i in [1..255]];assert num ne 0 and den ne 0;
out:="P79_P89/data/P79/quotient18/jmap_q_data.m";
PrintFile(out,Sprintf("N:=79;genus:=18;precision:=600;SturmBoundWeight16:=320;qcoeffs600:=Matrix(GF(2),18,600,%o);Diamond3Matrix:=Matrix(GF(2),18,18,%o);Std2Exponents:=%o;Std8Exponents:=%o;Gvec:=%o;NumeratorVec:=%o;DenominatorVec:=%o;",Eltseq(Q2),Eltseq(Md),[Exponents(m):m in S2],[Exponents(m):m in S8],Eltseq(gv),Eltseq(av),Eltseq(bv)):Overwrite:=true);
printf "NONIC79_JMAP_IDENTITIES numerator=E4cubed_times_G denominator=Delta_times_G precision=600 Sturm_bound=320 G_space_dimension=%o time=%o\n",Dimension(K),Cputime(tt);
PrintFile(out,"print \"NONIC79_QUOTIENT_JMAP_Q_DATA_READY\";");print "NONIC79_JMAP_Q_EXPORT_COMPLETE";quit;
