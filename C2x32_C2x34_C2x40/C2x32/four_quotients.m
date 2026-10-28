/* Complete necessary degree9 finite support filter, including all target
   cuspidal complements. Four rational maps to X1(32), no degree8 exclusion. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_diamond_hash_m2.m");
Attach("C2x32_C2x34_C2x40/C2x32/deps/fast_full_aut_hash_m2_32.m");
load "common/pair_af_place_linear.m";
Y:=MDX1(32,GF(3));C:=Curve(Y);assert Genus(C) eq 17;
os:=CuspOrbitsQ(Y);cs:=Cusps(Y);c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
F0,mp:=AlgorithmicFunctionField(FunctionField(C));F:=RationalExtensionRepresentation(F0);
base:=DivisorGroup(F)!FunctionFieldDivisor(Divisor(c0));G,_,phi:=ClassGroup(F:ReductionDivisor:=base,Proof:=true);
hg:=[phi(FunctionFieldDivisor(o-Degree(o)*Divisor(c0))):o in os];H:=sub<G|hg>;
assert Invariants(H) eq [2,2,2,4,120,11640];GG,pi:=quo<G|H>;
cc:=[phi(FunctionFieldDivisor(Divisor(c)-Degree(c)*Divisor(c0))):c in cs];
allowed:=AssociativeArray();for d in [0..4] do allowed[d]:={GG|};end for;allowed[0]:={GG!0};
for i in [1..#cs] do w:=Degree(cs[i]);for d in [w..4] do allowed[d] join:={z+pi(cc[i]):z in allowed[d-w]};end for;end for;
coords:=[F!mp(FunctionField(C).i):i in [1,2]];
print "C232_FOURQUOTIENT_GROUP",Invariants(G),Invariants(H),[#allowed[d]:d in [0..4]],Cputime(tt);
X:=MDX11(2,32,GF(3));assert Genus(X) eq 37;fp:=FrobeniusPolynomial(JH(64,[33]),3);traces:=[];closed:=[];
for n in [1..9] do
 sn:=-n*Coefficient(fp,74-n);for j in [1..n-1] do sn-:=Coefficient(fp,74-j)*traces[n-j];end for;Append(~traces,Integers()!sn);
 count:=3^n+1-Integers()!sn;for d in Divisors(n) do if d lt n then count-:=d*closed[d];end if;end for;assert count mod n eq 0;Append(~closed,count div n);
end for;
records:=[* *];coverage:=[* *];remaining:=[* *];alltarget:=[* *];cache:=AssociativeArray();R:=PolynomialRing(GF(3));
for d in [6..9] do
 xs:=NoncuspidalPlaces(X,d);assert #xs eq closed[d]-#[x:x in Cusps(X)|Degree(x) eq d];
 rs,labels,keys:=PlacesUpToFullAutM2ThirtyTwo(X,xs);K:=Parent(Evaluate(X`_coordinates[1],xs[1]));
 pol:=[Integers()!a:a in Coefficients(DefiningPolynomial(K,GF(3)))];
 ok:=[M2ThirtyTwoOrbitKeys(X,x,K):x in rs];Append(~coverage,[*d,pol,keys,labels,ok*]);
 print "C232_FOURQUOTIENT_TYPE_START",d,#xs,#rs,Cputime(tt);
 for i in [1..#rs] do
  bc:=MDM2ToCommonField([Evaluate(f,rs[i]):f in X`_coordinates],K);b,c:=Explode(bc);
  E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P:=E![0,0];Q:=E![b,b];assert Order(P) eq 2 and Order(Q) eq 32 and P ne 16*Q;
  maps:=[* *];passed:=true;
  for variant in [1..4] do
   if variant le 2 then EG:=E;U:=variant eq 1 select Q else P+Q;
   else
    T:=variant eq 3 select P else P+16*Q;EG,ph:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);U:=ph(Q);
   end if;assert Order(U) eq 32;
   pair:=TateNormalForm_xy(EG,U);key:=[d] cat MDM2FrobeniusKey(pair,K,3,d);
   if IsDefined(cache,key) then ti:=cache[key];cl:=alltarget[ti][3];e:=alltarget[ti][2];
   else
    yp:=C48PairAFPlaceLinear(coords,pair);e:=Degree(yp);assert d mod e eq 0;cl:=(d div e)*phi(Divisor(yp)-e*base);
    Append(~alltarget,[*key,e,cl*]);ti:=#alltarget;cache[key]:=ti;
   end if;
   flag:=-pi(cl) in allowed[9-d];Append(~maps,[*variant,ti,e,d div e,Eltseq(cl),flag*]);
   if not flag then passed:=false;end if;
  end for;
  vals:=[[Integers()!a:a in Eltseq(z,GF(3))]:z in bc];Append(~records,[*d,i,pol,vals,maps,passed*]);
  if passed then Append(~remaining,[*d,i,pol,vals*]);end if;
  print "C232_FOURQUOTIENT_CASE",d,i,passed,Cputime(tt);
 end for;
 Write("C2x32_C2x34_C2x40/data/C2x32/four_quotients_arithmetic.txt",Sprint([*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(g):g in hg],[Degree(c):c in cs],[Eltseq(a):a in cc],records,remaining,coverage,Coefficients(fp)*]):Overwrite:=true);
 Write("C2x32_C2x34_C2x40/data/C2x32/four_quotients_survivors.m","C232FourQuotientSurvivors:=" cat Sprint(remaining) cat ";":Overwrite:=true);
end for;
printf "C232_FOURQUOTIENT_DONE tested=%o remaining=%o complete=true cpu=%o\n",#records,#remaining,Cputime(tt);quit;
