SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");
assert assigned N and assigned RecordFile;N:=StringToInteger(N);
Ds:=[StringToInteger(s):s in Split(Ds,",")];Qs:=[StringToInteger(s):s in Split(Qs,",")];
if N in [41,47,49,59,71] then a:=1;else a:=PrimitiveRoot(N);end if;
if N eq 61 then a:=PrimitiveRoot(N)^6 mod N;end if;
t0:=Cputime();X:=MDX1(N,GF(2));C:=Curve(X);g:=Genus(C);if N eq 49 then assert g eq 69;else assert g eq (N-5)*(N-7) div 24;end if;
printf "NONIC_PRIME_MODEL N=%o genus=%o alpha_diamond=%o time=%o\n",N,g,a,Cputime(t0);
PrintFile(RecordFile,Sprintf("N:=%o;Ds:=%o;Qs:=%o;a:=%o;Inputs:=[**];",N,Ds,Qs,a):Overwrite:=true);
summary:=[];
for d in Ds do
 pl:=NoncuspidalPlaces(X,d);reps:=PlacesUpToDiamond(X,pl);
 printf "NONIC_PRIME_ENUM N=%o degree=%o all_places=%o diamond_reps=%o time=%o\n",N,d,#pl,#reps,Cputime(t0);
 survivors:=[];
 for i in [1..#reps] do
  P:=reps[i];D:=Divisor(P);E:=EllipticCurve(X,P);PN:=LevelStructure(X,P)`P;
  assert Order(PN) eq N;
  K:=BaseRing(E);assert Degree(K) eq d;
  PrintFile(RecordFile,Sprintf("Append(~Inputs,<%o,%o,%o,%o,%o>);",d,i,[Integers()!c:c in Eltseq(DefiningPolynomial(K,GF(2)))],[[Integers()!c:c in Eltseq(v,GF(2))]:v in aInvariants(E)],[[Integers()!c:c in Eltseq(PN[k],GF(2))]:k in [1..2]]));
  eliminated:=false;
  for q in Qs do
   assert IsPrime(q) and IsOdd(q) and N mod q ne 0;
   T:=HeckeOperatorFast(X,q,P);assert Degree(T) eq (q+1)*d;
   B:=T-D-q*Divisor(DiamondOperator(X,q,P));
   if a ne 1 then B:=DiamondOperator(X,a,B)-B;end if;
   assert Degree(B) eq 0;
   printf "NONIC_PRIME_TEST_START N=%o degree=%o index=%o q=%o time=%o\n",N,d,i,q,Cputime(t0);
   flag:=IsPrincipal(B);
   printf "NONIC_PRIME_TEST_RESULT N=%o degree=%o index=%o q=%o principal=%o time=%o\n",N,d,i,q,flag,Cputime(t0);
   if not flag then eliminated:=true;break;end if;
  end for;
  if not eliminated then Append(~survivors,i);end if;
 end for;
 printf "NONIC_PRIME_DEGREE_SUMMARY N=%o degree=%o all_places=%o diamond_reps=%o surviving_indices=%o time=%o\n",N,d,#pl,#reps,survivors,Cputime(t0);
 Append(~summary,<d,#pl,#reps,#survivors>);
end for;
printf "NONIC_PRIME_GOOD_BRANCHES_COMPLETED N=%o summary=%o time=%o\n",N,summary,Cputime(t0);quit;
