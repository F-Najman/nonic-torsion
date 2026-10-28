SetColumns(0);SetQuitOnError(true);SetMemoryLimit(5*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");Attach("common/fast_hecke.m");load "common/pair_af_place_linear.m";load "C48/data/hecke13_1/inputs.m";
X:=MDX1(48,GF(5));FC:=FunctionField(Curve(X));F,am:=AlgorithmicFunctionField(FC);coords:=[am(FC.i):i in [1..2]];assert Genus(F) eq 37 and 13 mod 4 eq 1;
fields:=AssociativeArray();for d in [3..8] do fields[d]:=ext<GF(5)|Polynomial(GF(5),FieldPolynomials[d-2])>;end for;
cache:=AssociativeArray();out:=Open("C48/data/hecke13_1/results.py","w");fprintf out,"q=13\ntests=[\n";remaining:=[];
for k in [1..#Cases] do
 degrees,indices,pairs:=Explode(Cases[k]);assert &+degrees eq 8;AA:=DivisorGroup(F)!0;
 for h in [1..#degrees] do
  d:=degrees[h];raw:=pairs[h];ky:=[d] cat &cat [[Integers()!v:v in part]:part in raw];
  if not IsDefined(cache,ky) then
   K:=fields[d];bc:=[Evaluate(PolynomialRing(GF(5))!v,K.1):v in raw];b,c:=Explode(bc);E:=EllipticCurve([1-c,-b,-b,K!0,K!0]);P:=E![0,0];assert 48*P eq E!0 and 24*P ne E!0 and 16*P ne E!0;
   pp:=C48PairAFPlaceLinear(coords,TateNormalForm_xy(E,P));assert Degree(pp) eq d;DD:=DivisorGroup(F)!0;
   images:=MDIsogeniesFast(E,P,13);printf "NONIC48_HEC13_INPUT case=%o degree=%o kernel_orbits=%o time=%o\n",CaseIDs[k],d,#images,Cputime(tt);
   for z in images do EG,PG,e:=Explode(z);image:=C48PairAFPlaceLinear(coords,TateNormalForm_xy(EG,PG));m:=e*d/Degree(image);assert Denominator(m) eq 1;DD+:=(Integers()!m)*Divisor(image);printf "NONIC48_HEC13_IMAGE case=%o source_degree=%o kernel_orbit=%o image_degree=%o weight=%o time=%o\n",CaseIDs[k],d,e,Degree(image),m,Cputime(tt);end for;
   assert Degree(DD) eq 14*d;dp:=C48PairAFPlaceLinear(coords,TateNormalForm_xy(E,13*P));assert Degree(dp) eq d;cache[ky]:=DD-Divisor(pp)-13*Divisor(dp);assert Degree(cache[ky]) eq 0;
  end if;AA+:=cache[ky];
 end for;
 printf "NONIC48_HEC13_RR_START case=%o time=%o\n",CaseIDs[k],Cputime(tt);bs:=Basis(AA:Reduction:=false,Simplification:="None");assert #bs in [0,1];if #bs ne 0 then Append(~remaining,CaseIDs[k]);end if;
 fprintf out,"[%o,%o,%o,%o],\n",CaseIDs[k],degrees,indices,#bs;Flush(out);printf "NONIC48_HEC13_CASE_DONE case=%o h0=%o time=%o\n",CaseIDs[k],#bs,Cputime(tt);
end for;
fprintf out,"]\nremaining=%o\ncomplete=True\n",remaining;delete out;
printf "NONIC48_HEC13_BATCH_DONE tested=%o remaining=%o direct_RR=true complete=true time=%o\n",#Cases,#remaining,Cputime(tt);quit;
