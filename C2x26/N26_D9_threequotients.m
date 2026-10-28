/* Complete NEW degree9 pure-support quotient filter on X1(2,26).
   Low-degree noncuspidal supports5..8 require separate cusp complement
   filters. This stage never discards a source solely by finite point count. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_diamond_hash_m2.m");
load "common/pair_af_place_linear.m";
X:=MDX11(2,26,GF(3));assert Genus(X) eq 25;
xs:=NoncuspidalPlaces(X,9);assert #xs eq 2328;
K:=Parent(Evaluate(X`_coordinates[1],xs[1]));assert #K eq 3^9;
seen:=AssociativeArray();reps:=[];inputkeys:=[];labels:=[];orbitkeys:=[* *];
for x in xs do
 bc:=MDM2CanonicalPair(X,x,K);key:=MDM2FrobeniusKey(bc,K,3,9);Append(~inputkeys,key);
 if IsDefined(seen,key) then Append(~labels,seen[key]);continue;end if;
 Append(~reps,x);Append(~labels,#reps);b,c:=Explode(bc);
 E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P:=E![0,0];Q:=E![b,b];
 twos:=[P,13*Q,P+13*Q];V:=14*Q;keys:=[];
 for ri in [1..3] do for pi in [1..3] do
  if ri eq pi then continue;end if;
  for a in [1..6] do
   U:=twos[pi];W:=twos[ri]+a*V;assert Order(U) eq 2 and Order(W) eq 26 and U ne 13*W;
   r:=U[1];xx:=W[1]-r;yy:=W[2];assert xx ne 0 and yy ne 0;
   bn:=xx^3/yy^2;cn:=(c+3*r)*xx^2/yy^2;
   ky:=MDM2FrobeniusKey([bn,cn],K,3,9);Append(~keys,ky);
   if IsDefined(seen,ky) then assert seen[ky] eq #reps;else seen[ky]:=#reps;end if;
  end for;
 end for;end for;assert #keys eq 36;Append(~orbitkeys,keys);
end for;
assert #Keys(seen) eq #xs;
print "N26_D9_THREEQUOTIENT_ORBITS",#xs,#reps,Cputime(tt);
Y:=MDX1(26,GF(3));assert Genus(Y) eq 10;CY:=Curve(Y);
os:=CuspOrbitsQ(Y);c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
F0,mp:=AlgorithmicFunctionField(FunctionField(CY));F:=RationalExtensionRepresentation(F0);
c0D:=F!1; // divisor below avoids generic curve-place conversion later
base:=DivisorGroup(F)!FunctionFieldDivisor(Divisor(c0));
G,_,phi:=ClassGroup(F:ReductionDivisor:=base,Proof:=true);
hg:=[phi(FunctionFieldDivisor(o-Degree(o)*Divisor(c0))):o in os];H:=sub<G|hg>;
assert Invariants(H) eq [133,1995];
fp:=FrobeniusPolynomial(JOne(26),3);assert &*[v:v in Invariants(G)|v ne 0] eq Integers()!Evaluate(fp,1);
coords:=[F!mp(FunctionField(CY).i):i in [1,2]];
print "N26_D9_THREEQUOTIENT_GROUP",Invariants(G),Invariants(H),Cputime(tt);
records:=[* *];remaining:=[* *];targetrecords:=[* *];targetkeys:=AssociativeArray();
for i in [1..#reps] do
 x:=reps[i];b,c:=Explode(MDM2CanonicalPair(X,x,K));
 E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P:=E![0,0];Q:=E![b,b];
 maps:=[* *];survives:=true;
 for variant in [1..3] do
  U:=[Q,P+Q,P+14*Q][variant];assert Order(U) eq 26;
  pair:=TateNormalForm_xy(E,U);key:=MDM2FrobeniusKey(pair,K,3,9);
  if IsDefined(targetkeys,key) then ti:=targetkeys[key];cl:=targetrecords[ti][3];e:=targetrecords[ti][2];
  else
   yp:=C48PairAFPlaceLinear(coords,pair);e:=Degree(yp);assert 9 mod e eq 0;
   cl:=(9 div e)*phi(Divisor(yp)-e*base);Append(~targetrecords,[*key,e,cl*]);ti:=#targetrecords;targetkeys[key]:=ti;
  end if;
  passed:=cl in H;Append(~maps,[*variant,ti,e,9 div e,Eltseq(cl),passed*]);if not passed then survives:=false;end if;
 end for;
 pol:=[Integers()!z:z in Coefficients(DefiningPolynomial(K,GF(3)))];vals:=[[Integers()!z:z in Eltseq(v,GF(3))]:v in [b,c]];
 Append(~records,[*i,pol,vals,maps,survives*]);
 if survives then Append(~remaining,[*i,pol,vals*]);end if;
 print "N26_D9_THREEQUOTIENT_CASE",i,survives,Cputime(tt);
end for;
Write("C2x26/data/N26_D9_threequotients_arithmetic.txt",Sprint([*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(g):g in hg],records,remaining,Coefficients(fp),inputkeys,labels,orbitkeys*]):Overwrite:=true);
Write("C2x26/data/N26_D9_threequotients_survivors.m","N26D9ThreeQuotientSurvivors:=" cat Sprint(remaining) cat ";":Overwrite:=true);
printf "N26_D9_THREEQUOTIENT_DONE tested=%o remaining=%o complete=true time=%o\n",#records,#remaining,Cputime(tt);quit;
