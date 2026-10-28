SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_diamond_hash_m2.m");
Attach("C2x28/deps/fast_full_aut_hash_m2_28.m");
load "C2x28/data/four_quotients_survivors.m";
R:=PolynomialRing(GF(3));output:=[* *];
for rec in C228FourQuotientSurvivors do
 d,i,pol,vals:=Explode(rec);K<z>:=ext<GF(3)|R!pol>;bc:=[Evaluate(R!v,z):v in vals];b,c:=Explode(bc);
 E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P:=E![0,0];Q:=E![b,b];key:=MDM2FrobeniusKey(bc,K,3,d);stabs:=[* *];
 for par in M2TwentyEightAutParameters() do
  b1,c1,a:=Explode(par);U:=P+14*b1*Q;V:=c1*P+a*Q;pair:=M2TwentyEightNormalizedPair(E,U,V);
  if MDM2FrobeniusKey(pair,K,3,d) eq key then
   phases:=[k:k in [0..d-1]|pair eq [v^(3^k):v in bc]];assert #phases eq 1;Append(~stabs,par cat [phases[1]]);
  end if;
 end for;
 Append(~output,[*d,i,stabs*]);print "C228_EXACT_STABILIZER",d,i,stabs;
end for;
Write("C2x28/data/stabilizers_arithmetic.txt",Sprint(output):Overwrite:=true);
print "C228_STABILIZERS_DONE complete=true";quit;
