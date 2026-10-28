// Root-allocated first canonical/unit stage only; no function-field normalization.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
N:=65;g:=121;QPrec:=1100;TatePrec:=1500;
MS:=CuspidalSubspace(ModularSymbolsH(N,[1],2,1));assert Dimension(MS)eq g;decomp:=NewformDecomposition(MS);
printf "NONIC65_ALLBAD_CANONICAL_UNIT_NEWFORMS dimensions=%o time=%o\n",[Dimension(A):A in decomp],Cputime(tt);
function Flat(z)
 if Type(z)eq FldRatElt or Type(z)eq RngIntElt then return [Rationals()|z];end if;return &cat[Flat(y):y in Eltseq(z)];
end function;
Q0:=[];D0:=[];opns:=[2,3,7,11,17];OP0:=AssociativeArray();aq:=AssociativeArray();for nn in opns do OP0[nn]:=[];end for;
for A in decomp do
 ff:=qEigenform(A,QPrec);KF:=BaseRing(Parent(ff));chi:=(Coefficient(ff,3)^2-Coefficient(ff,9))/3;
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
  if LRatio(A,1)eq 0 then assert lev eq 65 and Dimension(A)eq 1;for nn in opns do aq[nn]:=Rationals()!Coefficient(ff,nn);end for;end if;
 end if;
 printf "NONIC65_ALLBAD_CANONICAL_UNIT_EIGENFORM dimension=%o accumulated=%o time=%o\n",Dimension(A),#Q0,Cputime(tt);
end for;
QM0:=Matrix(Rationals(),Q0);DM0:=Matrix(Rationals(),D0);assert Nrows(QM0)eq g;
den:=LCM([Denominator(a):a in Eltseq(QM0)]);QMZ:=Saturation(ChangeRing(den*QM0,Integers()));QMR:=ChangeRing(QMZ,Rationals());assert Nrows(QMZ)eq g and RowSpace(QMR)eq RowSpace(QM0);
Q2:=ChangeRing(QMZ,GF(2));assert Rank(Q2)eq g and Rank(Submatrix(Q2,1,1,g,337))eq g;
TT:=Solution(QM0,QMR);assert TT*QM0 eq QMR;Dfull:=TT*DM0;MdQ:=Solution(QMR,Dfull);assert MdQ*QMR eq Dfull and &and[IsOdd(Denominator(a)):a in Eltseq(MdQ)];Md:=ChangeRing(MdQ,GF(2));assert Md^12 eq IdentityMatrix(GF(2),g)and Md ne IdentityMatrix(GF(2),g);
printf "NONIC65_ALLBAD_CANONICAL_UNIT_SATURATED_BASIS rank121=true precision1100=true Sturm_weight2=336 diamond3order12=true time=%o\n",Cputime(tt);
load "C65/data/source_units/canonical_unit_data.m";Ualign:=Solution(Q2,qcoeffs);assert Rank(Ualign)eq g and Ualign*Q2 eq qcoeffs and Ualign*Md*Ualign^-1 eq Diamond2Matrix;
OP:=AssociativeArray();
for nn in opns do OO:=Matrix(Rationals(),OP0[nn]);MM:=Solution(QMR,TT*OO);assert MM*QMR eq TT*OO and &and[IsOdd(Denominator(v)):v in Eltseq(MM)];OP[nn]:=Ualign*ChangeRing(MM,GF(2))*Ualign^-1;end for;
II:=IdentityMatrix(GF(2),g);allops:=[];meta:=[];
for nn in opns,ell in[3,7,11,17]do
 kk:=[<aa,bb>:aa in[0..1],bb in[0..11]|(2^aa*3^bb)mod 65 in[ell mod 65,(-ell)mod 65]];assert #kk eq 1;delta:=Diamond1Matrix^kk[1][1]*Diamond2Matrix^kk[1][2];AA:=(OP[nn]-(GF(2)!aq[nn])*II)*(OP[ell]-II-delta);Append(~allops,AA);Append(~meta,<nn,ell,aq[nn],kk[1],Rank(AA)>);
end for;
WW:=BasisMatrix(RowSpace(VerticalJoin(allops)));assert Nrows(WW)le 120;
PrintFile("C65/data/allbad/annihilator_data.m",Sprintf("N:=65;genus:=121;HeckeIndices:=%o;PositiveFactorEigenvalues:=%o;HeckeMatrices:=%o;AnnihilatorMetadata:=%o;AnnihilatingDifferentials:=Matrix(GF(2),%o,121,%o);",opns,[aq[n]:n in opns],[Eltseq(OP[n]):n in opns],meta,Nrows(WW),Eltseq(WW)):Overwrite:=true);
printf "NONIC65_ALLBAD_CANONICAL_ANNIHILATOR_COMPLETE dimensions=%o positive_eigenvalues=%o meta=%o time=%o\n",Nrows(WW),[aq[n]:n in opns],meta,Cputime(tt);quit;
