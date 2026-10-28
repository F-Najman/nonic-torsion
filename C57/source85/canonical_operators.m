// Root-allocated first canonical/unit stage only; no function-field normalization.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
N:=57;g:=85;QPrec:=1100;TatePrec:=1500;
MS:=CuspidalSubspace(ModularSymbolsH(N,[1],2,1));assert Dimension(MS)eq g;decomp:=NewformDecomposition(MS);
printf "NONIC57_SOURCE85_CANONICAL_UNIT_NEWFORMS dimensions=%o time=%o\n",[Dimension(A):A in decomp],Cputime(tt);
function Flat(z)
 if Type(z)eq FldRatElt or Type(z)eq RngIntElt then return [Rationals()|z];end if;return &cat[Flat(y):y in Eltseq(z)];
end function;
Q0:=[];D0:=[];opns:=[2,5,7,11,13];OP0:=AssociativeArray();aq:=AssociativeArray();for nn in opns do OP0[nn]:=[];end for;
for A in decomp do
 ff:=qEigenform(A,QPrec);KF:=BaseRing(Parent(ff));chi:=(Coefficient(ff,5)^2-Coefficient(ff,25))/5;
 block:=[];dblock:=[];
 if KF cmpeq Rationals()then
  Append(~block,[Rationals()!Coefficient(ff,n):n in[0..QPrec-1]]);Append(~dblock,[Rationals()!(chi*Coefficient(ff,n)):n in[0..QPrec-1]]);
 else
  dd:=#Flat(KF!1);cofs:=[Flat(KF!Coefficient(ff,n)):n in[0..QPrec-1]];dcofs:=[Flat(KF!(chi*Coefficient(ff,n))):n in[0..QPrec-1]];
  assert &and[#c eq dd:c in cofs]and &and[#c eq dd:c in dcofs];for k in[1..dd]do Append(~block,[c[k]:c in cofs]);Append(~dblock,[c[k]:c in dcofs]);end for;
 end if;
 lev:=Level(AssociatedNewSpace(A));assert N mod lev eq 0;
 for oldd in Divisors(N div lev)do
  for j in [1..#block]do
   Append(~Q0,[n mod oldd eq 0 select block[j][n div oldd+1]else Rationals()!0:n in[0..QPrec-1]]);
   Append(~D0,[n mod oldd eq 0 select dblock[j][n div oldd+1]else Rationals()!0:n in[0..QPrec-1]]);
  end for;
 end for;
 for nn in opns do
  cf:=[Flat(KF!(Coefficient(ff,nn)*Coefficient(ff,n))):n in[0..QPrec-1]];
  for oldd in Divisors(N div lev),k in[1..#block]do Append(~OP0[nn],[n mod oldd eq 0 select cf[n div oldd+1][k]else Rationals()!0:n in[0..QPrec-1]]);end for;
 end for;
 if KF cmpeq Rationals() then
  if LRatio(A,1)eq 0 then assert lev eq 57 and Dimension(A)eq 1;for nn in opns do aq[nn]:=Rationals()!Coefficient(ff,nn);end for;end if;
 end if;
 printf "NONIC57_SOURCE85_CANONICAL_UNIT_EIGENFORM dimension=%o accumulated=%o time=%o\n",Dimension(A),#Q0,Cputime(tt);
end for;
QM0:=Matrix(Rationals(),Q0);DM0:=Matrix(Rationals(),D0);assert Nrows(QM0)eq g;
den:=LCM([Denominator(a):a in Eltseq(QM0)]);QMZ:=Saturation(ChangeRing(den*QM0,Integers()));QMR:=ChangeRing(QMZ,Rationals());assert Nrows(QMZ)eq g and RowSpace(QMR)eq RowSpace(QM0);
Q2:=ChangeRing(QMZ,GF(2));assert Rank(Q2)eq g and Rank(Submatrix(Q2,1,1,g,241))eq g;
TT:=Solution(QM0,QMR);assert TT*QM0 eq QMR;Dfull:=TT*DM0;MdQ:=Solution(QMR,Dfull);assert MdQ*QMR eq Dfull and &and[IsOdd(Denominator(a)):a in Eltseq(MdQ)];Md:=ChangeRing(MdQ,GF(2));assert Md^18 eq IdentityMatrix(GF(2),g)and Md ne IdentityMatrix(GF(2),g);
printf "NONIC57_SOURCE85_CANONICAL_UNIT_SATURATED_BASIS rank85=true precision1100=true Sturm_weight2=240 diamond5order18=true time=%o\n",Cputime(tt);
load "C57/data/source85/qbasis_data.m";assert qcoeffs eq Q2 and Diamond5Matrix eq Md;
OP:=AssociativeArray();
for nn in opns do OO:=Matrix(Rationals(),OP0[nn]);MM:=Solution(QMR,TT*OO);assert MM*QMR eq TT*OO and &and[IsOdd(Denominator(v)):v in Eltseq(MM)];OP[nn]:=ChangeRing(MM,GF(2));end for;
II:=IdentityMatrix(GF(2),g);allops:=[];meta:=[];
for nn in opns,ell in[5,7,11,13]do
 kk:=[a:a in[0..17]|5^a mod 57 in[ell mod 57,(-ell)mod 57]];assert #kk eq 1;delta:=Md^kk[1];AA:=(OP[nn]-(GF(2)!aq[nn])*II)*(OP[ell]-II-delta);Append(~allops,AA);Append(~meta,<nn,ell,aq[nn],kk[1],Rank(AA)>);
end for;
WW:=BasisMatrix(RowSpace(VerticalJoin(allops)));assert Nrows(WW)le 84;
PrintFile("C57/data/source85/annihilator_data.m",Sprintf("N:=57;genus:=85;HeckeIndices:=%o;PositiveFactorEigenvalues:=%o;HeckeMatrices:=%o;AnnihilatorMetadata:=%o;AnnihilatingDifferentials:=Matrix(GF(2),%o,85,%o);",opns,[aq[n]:n in opns],[Eltseq(OP[n]):n in opns],meta,Nrows(WW),Eltseq(WW)):Overwrite:=true);
printf "NONIC57_SOURCE85_CANONICAL_ANNIHILATOR_COMPLETE dimensions=%o positive_eigenvalues=%o meta=%o time=%o\n",Nrows(WW),[aq[n]:n in opns],meta,Cputime(tt);quit;
