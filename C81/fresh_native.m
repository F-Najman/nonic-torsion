SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);
t0:=Cputime(); N:=81;
D:=FullDirichletGroup(N);total:=0;
for chi in GaloisConjugacyRepresentatives(D) do
 if IsOdd(chi) then continue;end if;
 chimin:=MinimalBaseRingCharacter(chi);
 S:=CuspidalSubspace(ModularSymbols(chimin,2,1));
 d:=Dimension(S); w:=Dimension(WindingSubmodule(S));
 assert d eq w;
 total+:=d*Degree(BaseRing(S));
 printf "NONIC81_WINDING order=%o conductor=%o dimension=%o winding=%o coefficient_degree=%o time=%o\n",Order(chi),Conductor(chi),d,w,Degree(BaseRing(S)),Cputime(t0);
end for;
assert total eq 190;
print "NONIC81_RANK_ZERO",total;
J:=JOne(81);assert Sign(J) eq 0 and Dimension(J) eq 190;
H:=IntegralHomology(J);assert Rank(H) eq 380;
d:=Rank(H);K:=VectorSpace(GF(2),d);
for q in [5,7] do
 T:=ChangeRing(IntegralMatrix(HeckeOperator(J,q)),Integers());
 Tqq:=ChangeRing(IntegralMatrix(HeckeOperator(J,q^2)),Integers());
 Dq:=ChangeRing((ChangeRing(T,Rationals())^2-ChangeRing(Tqq,Rationals()))/q,Integers());
 A:=T-IdentityMatrix(Integers(),d)-q*Dq;
 K:=K meet Nullspace(ChangeRing(A,GF(2)));
 printf "NONIC81_RATIONAL_TWO q=%o integral_homology_rank=%o joint_kernel_dimension=%o time=%o\n",q,d,Dimension(K),Cputime(t0);
end for;
assert Dimension(K) eq 0;
print "NONIC81_NATIVE_DONE";
quit;
