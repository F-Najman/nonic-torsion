SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
R<x>:=PolynomialRing(Rationals());PP:=R!1;total:=0;
for chi0 in GaloisConjugacyRepresentatives(FullDirichletGroup(133)) do
 if IsOdd(chi0) or Evaluate(chi0,78) ne 1 then continue;end if;
 chi:=MinimalBaseRingCharacter(chi0);S:=CuspidalSubspace(ModularSymbols(chi,2,1));W:=WindingSubmodule(S);K:=BaseRing(S);d:=Dimension(S);w:=Dimension(W);total+:=d*Degree(K);
 printf "Z133_WINDING order=%o conductor=%o dimension=%o winding=%o field_degree=%o\n",Order(chi),Conductor(chi),d,w,Degree(K);
 if d eq w then continue;end if;
 BS:=Matrix([Eltseq(b):b in Basis(S)]);BW:=Matrix([Eltseq(b):b in Basis(W)]);T:=HeckeOperator(S,2);assert T^2-HeckeOperator(S,4) eq 2*Evaluate(chi,2)*IdentityMatrix(K,d);
 printf "Z133_DIAMOND_RECURRENCE order=%o scalar=%o exact=true\n",Order(chi),Evaluate(chi,2);VS:=VectorSpace(K,d);
 ww:=[VS!Solution(BS,r):r in Rows(BW)];UW:=sub<VS|ww>;assert &and[r*T in UW:r in ww];bb:=[VS!v:v in Basis(UW)];
 for v in Basis(VS) do if v notin sub<VS|bb> then Append(~bb,v);end if;end for;
 B:=Matrix(bb);M:=B*T*B^-1;QQ:=Submatrix(M,w+1,w+1,d-w,d-w);degK:=Degree(K);dim:=(d-w)*degK;
 if Type(K) eq FldRat then MQ:=ChangeRing(QQ,Rationals());
 else
  rb:=Basis(K);rr:=[];
  for i in [1..d-w],b in rb do Append(~rr,&cat[Eltseq(b*QQ[i,j]):j in [1..d-w]]);end for;
  MQ:=Matrix(Rationals(),rr);
 end if;
 assert Nrows(MQ) eq dim and Ncols(MQ) eq dim;
 if Order(chi) eq 1 then
 assert Evaluate(x^4+4*x^3+x^2-8*x-3,MQ) eq 0;
 printf "Z133_TRIVIAL_FINITE_ANNIHILATOR coefficients=[-3,-8,1,4,1] exact=true\n";
end if;
if Order(chi) eq 3 then
 printf "Z133_POSITIVE_LINEAR_CHECK a2=%o chi2=%o ratio=%o\n",QQ[1,1],Evaluate(chi,2),QQ[1,1]/2;
 assert QQ[1,1] eq 2*Evaluate(chi,2)^2;
end if;
f:=R!MinimalPolynomial(MQ^3);PP:=LCM(PP,f);
 printf "Z133_POSITIVE_T2 order=%o rational_dimension=%o matrix=%o cubic_minimal_polynomial=%o\n",Order(chi),dim,Eltseq(MQ),Coefficients(f);
end for;
assert total eq 35;assert &and[Denominator(c) eq 1:c in Coefficients(PP)];
printf "Z133_SPLIT_FINITE_OPERATORS_VERIFIED coefficients=%o degree=%o genus=35 time=%o complete=true\n",Coefficients(PP),Degree(PP),Cputime(tt);
quit;
