SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
for chi0 in GaloisConjugacyRepresentatives(FullDirichletGroup(247)) do
 if IsOdd(chi0) or Evaluate(chi0,40) ne 1 or Order(chi0) eq 1 then continue;end if;
 chi:=MinimalBaseRingCharacter(chi0);q:=Order(chi) eq 2 select 17 else 3;
 S:=CuspidalSubspace(ModularSymbols(chi,2,1));K:=BaseRing(S);d:=Dimension(S);dk:=Degree(K);
 T:=HeckeOperator(S,q);E:=(1+q*Evaluate(chi,q))*IdentityMatrix(K,d)-T;z:=Determinant(E);
 if Type(K) eq FldRat then MQ:=ChangeRing(E,Rationals());n:=Integers()!z;fp:=[0,1];
 else
  rr:=[];for i in [1..d],b in Basis(K) do Append(~rr,&cat[Eltseq(b*E[i,j]):j in [1..d]]);end for;
  MQ:=Matrix(Rationals(),rr);n:=Integers()!Norm(z);fp:=Coefficients(DefiningPolynomial(K));
 end if;
 assert Nrows(MQ) eq d*dk and Ncols(MQ) eq d*dk;
 assert Determinant(MQ) eq n and IsOdd(n);
 file:=Sprintf("C247/data/norm_order%o_q%o.json",Order(chi),q);out:=Open(file,"w");
 fprintf out,"{\"order\":%o,\"q\":%o,\"relative_dimension\":%o,\"coefficient_degree\":%o,\"field_polynomial\":%o,\"rational_matrix_dimension\":%o,\"entries\":%o,\"exact_determinant_norm\":%o,\"odd\":true}\n",Order(chi),q,d,dk,fp,Nrows(MQ),[[Numerator(a),Denominator(a)]:a in Eltseq(MQ)],n;
 delete out;
 printf "NONIC247_NORM_EXPORTED order=%o q=%o dimension=%o norm=%o odd=true time=%o\n",Order(chi),q,d*dk,n,Cputime(tt);
end for;
print "NONIC247_NORM_EXPORT_DONE";quit;
