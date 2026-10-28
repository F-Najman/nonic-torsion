SetColumns(0);SetQuitOnError(true);SetMemoryLimit(10^9);tt:=Cputime();AttachSpec("mdmagma/v2/mdmagma.spec");load "C48/degree8/two_isogeny_48_to_2_24.m";load "C48/data/parameter9/records.m";load "C48/data/source9_survivors/inputs.m";
R:=PolynomialRing(GF(5));K<z>:=ext<GF(5)|R!FieldPolynomial>;function Key(bc)return Minimum([Eltseq(bc[1]^(5^i),GF(5)) cat Eltseq(bc[2]^(5^i),GF(5)):i in [0..8]]);end function;
function AutKey(E,P2,Q24)
 return Minimum([Key(DerickxNormalForm_bc(E,P2+eps*(12*Q24),a*Q24+eta*P2)):eps,eta in [0,1],a in [1,5,7,11]]);
end function;
seen:={};records:=[* *];expectedWeight:=0;
for idx in PassingIndices do
 rec:=Records[idx];bc:=[Evaluate(R!rec[j],z):j in [2,3]];b,c:=Explode(bc);E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P2:=E![0,0];Q24:=E![b,b];localTarget:={};localSource:={};expectedWeight+:=rec[5];
 for eps,eta in [0,1],a in [1,5,7,11] do
  pair:=DerickxNormalForm_bc(E,P2+eps*(12*Q24),a*Q24+eta*P2);key:=Key(pair);if key in localTarget then continue;end if;Include(~localTarget,key);bb,cc:=Explode(pair);yes,root:=IsSquare(bb);if not yes then continue;end if;
  dd:=bb*(1-bb-cc);AA:=-2*cc;BB:=cc^2-4*dd;ES:=EllipticCurve([K!0,AA,K!0,BB,K!0]);
  for sign in [K!1,-K!1] do
   xx:=cc+2*bb+2*sign*root;yy:=8*bb*xx^2/(BB-xx^2);PP:=ES![xx,yy];assert 48*PP eq ES!0 and 24*PP ne ES!0 and 16*PP ne ES!0;sk:=Key(TateNormalForm_bc(ES,PP));assert sk notin localSource;Include(~localSource,sk);
   orbit:={Key(TateNormalForm_bc(ES,u*PP)):u in [1,5,7,11,13,17,19,23]};assert #orbit eq 8;canonical:=Minimum(orbit);if canonical in seen then continue;end if;Include(~seen,canonical);
   cb,ccc:=Explode([Evaluate(R!canonical[j..j+8],z):j in [1,10]]);EC:=EllipticCurve([1-ccc,-cb,-cb,K!0,K!0]);PC:=EC![0,0];assert 48*PC eq EC!0 and 24*PC ne EC!0 and 16*PC ne EC!0;ET,LC:=C48TwoIsogenousModuli(EC,PC);assert AutKey(ET,LC`P,LC`Q) eq rec[2] cat rec[3];
   Append(~records,[*idx,canonical[1..9],canonical[10..18],8*]);
  end for;
 end for;assert #localSource eq rec[5];
end for;
assert #records*8 eq expectedWeight;out:=Open("C48/data/source9_survivors/records.m","w");fprintf out,"Source9FieldPolynomial:=%o;\nSource9Cases:=%o;\nComplete:=true;\n",FieldPolynomial,records;delete out;
printf "NONIC48_SOURCE9_EXPANSION_DONE targetorbits=%o sourceDiamondOrbits=%o sourceclosed=%o disjoint=true exact48=true complete=true time=%o\n",#PassingIndices,#records,expectedWeight,Cputime(tt);quit;
