SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);
H:={1,7,11};assert {a*b mod 19:a,b in H} eq H;
for e in [1..3] do
 F:=GF(2^e);seen:=0;
 for jj in F do
  for E in Twists(EllipticCurveFromjInvariant(jj)) do
   if #E mod 7 ne 0 then continue;end if;
   seen+:=1;t:=2^e+1-#E;vals:=[(a*a-t*a+2^e) mod 19:a in Sort(Setseq(H))];
   assert 0 notin vals;
   printf "NONIC133_Y_CUBIC_LOCAL e=%o Eorder=%o trace=%o H19=%o characteristic_values=%o\n",e,#E,t,Sort(Setseq(H)),vals;
  end for;
 end for;
 printf "NONIC133_Y_NO_GOOD e=%o qualifying_P7_curves=%o exact_all_j_and_twists=true\n",e,seen;
end for;
print "NONIC133_Y_ALLBAD_CUBIC_PASS";quit;
