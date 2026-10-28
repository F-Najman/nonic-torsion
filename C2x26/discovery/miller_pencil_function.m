/* Explicit square root by elliptic line identities, no root extraction. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
if not assigned Primes then Primes:=[3,7];end if;
for prime in Primes do
 kk:=prime eq 0 select Rationals() else GF(prime);X:=MDX11(2,26,kk);C:=Curve(X);FF:=FunctionField(C);
 print "C226_MILLER_MODEL_READY",prime,Cputime(tt);
 E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];b:=Q[1];c:=aInvariants(E)[2];dd:=aInvariants(E)[4];assert Q[2] eq b and P eq E![0,0];
 ts:=[P,13*Q,P+13*Q];print "C226_MILLER_TWO_TORSION_READY",prime,Cputime(tt);
 exps:=eval Read("C2x26/data/discovery/weighted_squared_unit_coefficients.txt");assert #exps eq 57 and exps[1..3] eq [0,0,0];
 assert &and[exps[4+4*(k-1)] eq 0:k in [1..12]];
 odd:=[&+[exps[4+4*(k-1)+j]:k in [1,3,5,7,9,11]]:j in [1..3]];assert odd eq [-1,1,1];
 result:=FF!1;points:=[k*Q:k in [1..8]];
 for j in [1..3] do
  r:=ts[j][1];a:=Q[1]-r;s:=FF!1;
  for k in [1..8] do
   e:=exps[4+4*(k-1)+j];if e ne 0 then result*:=s^e;end if;
   assert points[k][1]-r eq a^(k mod 2)*s^2;
   if k lt 8 then
    R:=points[k];if k eq 1 then slope:=(3*b^2+2*c*b+dd)/(2*b);else slope:=(R[2]-Q[2])/(R[1]-Q[1]);end if;
    line:=R[2]+slope*(r-R[1]);s:=line/(s*a);
   end if;
  end for;
 end for;
 planes:=[FF.1,FF.2,FF.1+1,FF.2+1,FF.1-FF.2,FF.1+FF.2+2];for i in [1..6] do e:=exps[51+i];assert e mod 2 eq 0;if e ne 0 then result*:=planes[i]^(e div 2);end if;end for;
 print "C226_MILLER_EXPLICIT_FUNCTION_READY",prime,Cputime(tt);
 Write("C2x26/data/discovery/miller_function_F" cat IntegerToString(prime) cat ".txt",Sprint(result):Overwrite:=true);
 if prime ne 0 then
  AF,am:=AlgorithmicFunctionField(FF);f:=am(result);cs:=[cu:cu in Cusps(X)|Degree(cu) eq 1];bv:=[Valuation(b,cu):cu in cs];
  Avals:=[3,-2,4,-18,1,6,-10,11,10];Amults:=[2,2,2,1,1,1,1,1,1];Cvals:=[2,5,-14,-6,9,12];Cmults:=[3,3,3,1,1,1];
  A:=&+[Amults[i]*Divisor(cs[Index(bv,Avals[i])]):i in [1..#Avals]];CC:=&+[Cmults[i]*Divisor(cs[Index(bv,Cvals[i])]):i in [1..#Cvals]];
  assert Divisor(f) eq FunctionFieldDivisor(A-CC);cu:=cs[Index(bv,7)];pp:=Support(FunctionFieldDivisor(Divisor(cu)))[1];value:=Evaluate(f,pp);print "C226_MILLER_PENCIL_DIVISOR_AND_VALUE_PASS",prime,value,Cputime(tt);
 end if;
end for;print "C226_MILLER_PENCIL_FUNCTION_DONE complete=true",Cputime(tt);quit;
