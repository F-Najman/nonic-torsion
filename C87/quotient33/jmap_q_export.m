SetColumns(0);SetQuitOnError(true);SetMemoryLimit(32*10^9);tt:=Cputime();
load "C87/data/quotient33/canonical_data.m";
gg:=33;prec:=900;MS:=CuspidalSubspace(ModularSymbolsH(87,[16],2,1));decomp:=NewformDecomposition(MS);
function Flat(z)
 if Type(z) eq FldRatElt or Type(z) eq RngIntElt then return [Rationals()|z];end if;
 return &cat[Flat(y):y in Eltseq(z)];
end function;
rows:=[];drows:=[];
for A in decomp do
 lev:=Level(AssociatedNewSpace(A));assert lev in [29,87];ff:=qEigenform(A,prec);KF:=BaseRing(Parent(ff));first:=#rows+1;chi:=(Coefficient(ff,2)^2-Coefficient(ff,4))/2;
 if KF cmpeq Rationals() then
  Append(~rows,[Rationals()!Coefficient(ff,n):n in [0..prec-1]]);Append(~drows,[Rationals()!(chi*Coefficient(ff,n)):n in [0..prec-1]]);
 else
  dd:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(ff,n)):n in [0..prec-1]];dcofs:=[Flat(KF!(chi*Coefficient(ff,n))):n in [0..prec-1]];
  for k in [1..dd] do Append(~rows,[z[k]:z in cofs]);Append(~drows,[z[k]:z in dcofs]);end for;
 end if;
 if lev eq 29 then
  oldRows:=rows[first..#rows];oldD:=drows[first..#drows];
  for r in oldRows do Append(~rows,[n mod 3 eq 0 select r[n div 3+1] else Rationals()!0:n in [0..prec-1]]);end for;
  for r in oldD do Append(~drows,[n mod 3 eq 0 select r[n div 3+1] else Rationals()!0:n in [0..prec-1]]);end for;
 end if;
end for;assert #rows eq gg and #drows eq gg;
Q0:=Matrix(Rationals(),rows);D0:=Matrix(Rationals(),drows);old:=ChangeRing(qcoeffs,Rationals());T:=Solution(Submatrix(Q0,1,1,gg,190),old);full:=T*Q0;dfull:=T*D0;assert Submatrix(full,1,1,gg,190) eq old;
assert &and[IsOdd(Denominator(a)):a in Eltseq(full)];Q2:=ChangeRing(full,GF(2));MdQ:=Solution(old,Submatrix(dfull,1,1,gg,190));assert MdQ*old eq Submatrix(dfull,1,1,gg,190);Md:=ChangeRing(MdQ,GF(2));assert Md^4 eq IdentityMatrix(GF(2),gg) and Md^2 ne IdentityMatrix(GF(2),gg);
printf "NONIC87_JMAP_Q_BASIS_READY precision900=true diamond2_order4=true time=%o\n",Cputime(tt);
P:=ProjectiveSpace(GF(2),gg-1);R:=CoordinateRing(P);mons:=[R.(z[1])*R.(z[2]):z in pairs];eqs:=[&+[v[i]*mons[i]:i in [1..#pairs]]:v in Rows(QuadRelations)];C:=Curve(P,eqs);IC:=Ideal(C);GB:=GroebnerBasis(IC);LM:=[LeadingMonomial(f):f in GB];
ls:=[&+[Md[i,j]*R.j:j in [1..gg]]:i in [1..gg]];assert &and[NormalForm(Evaluate(f,ls),IC) eq 0:f in eqs];
function StdMons(dd)
 ss:={R!1};for k in [1..dd] do cand:={m*R.i:m in ss,i in [1..33]};ss:={m:m in cand|not exists{l:l in LM|IsDivisibleBy(m,l)}};assert #ss eq (k eq 1 select 33 else (2*k-1)*32);end for;return Sort(Setseq(ss));
end function;
S3:=StdMons(3);S9:=StdMons(9);assert #S3 eq 160 and #S9 eq 544;
PS<q>:=PowerSeriesRing(GF(2),prec);Q:=[PS!Eltseq(v)+O(q^prec):v in Rows(Q2)];
function ExpMat(ms)return Matrix(GF(2),[[Coefficient(f,n):n in [0..899]] where f:=Evaluate(m,Q):m in ms]);end function;
E3:=ExpMat(S3);E9:=ExpMat(S9);assert Rank(E3) eq 160 and Rank(E9) eq 544;K:=Kernel(VerticalJoin(E3,E9));assert Dimension(K) eq 64;
v:=K.1;gv:=Vector(GF(2),Eltseq(v)[1..160]);av:=Vector(GF(2),Eltseq(v)[161..704]);assert gv ne 0 and gv*E3 eq av*E9;
G:=PS!Eltseq(gv*E3)+O(q^prec);Delta:=q*&*[(1+q^(8*n))*(1+q^(16*n)):n in [1..112]]+O(q^prec);bv:=Solution(E9,Vector(GF(2),[Coefficient(Delta*G,n):n in [0..899]]));assert bv*E9 eq Vector(GF(2),[Coefficient(Delta*G,n):n in [0..899]]);
PrintFile("C87/data/quotient33/jmap_q_data.m",Sprintf("N:=87;genus:=33;precision:=900;SturmBoundWeight18:=720;qcoeffs900:=Matrix(GF(2),33,900,%o);Diamond2Matrix:=Matrix(GF(2),33,33,%o);Std3Exponents:=%o;Std9Exponents:=%o;Gvec:=%o;NumeratorVec:=%o;DenominatorVec:=%o;",Eltseq(Q2),Eltseq(Md),[Exponents(m):m in S3],[Exponents(m):m in S9],Eltseq(gv),Eltseq(av),Eltseq(bv)):Overwrite:=true);
printf "NONIC87_JMAP_Q_CERTIFIED Gdegree3=true numerator_denominator_degree9=true weight18_Sturm720=true precision900=true Gspace_dimension64=true time=%o\n",Cputime(tt);quit;
