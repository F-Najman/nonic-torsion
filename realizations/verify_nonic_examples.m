/* Exact full-torsion certificates for all19 archived vanHoeij nonic places.
   No earlier-degree conclusion or dependency is used. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
load "common/data/inputs/vanhoeij_nonic_input.m";
out9:=Open("realizations/logs/verify_nonic_examples.log","w");
procedure Note9(s)
 fprintf out9,"%o\n",s;Flush(out9);print s;
end procedure;
R<T>:=PolynomialRing(Rationals());
function FromXY(xx,yy)
 r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx);
 s:=(xx*yy-yy+1)/(xx*yy);b:=r*s*(r-1);c:=s*(r-1);
 return EllipticCurve([1-c,-b,-b,0,0]);
end function;
function EvalMod(f,a)
 L:=Parent(a);
 return &+[L!Numerator(Coefficient(f,i))/(L!Denominator(Coefficient(f,i)))*a^i:i in [0..Degree(f)]];
end function;
for idx in [1..#NonicExamples] do
 N,fc,xc,yc:=Explode(NonicExamples[idx]);
 f:=R!fc;xp:=R!xc;yp:=R!yc;
 assert Degree(f) eq 9 and IsIrreducible(f) and IsMonic(f);
 K<a>:=NumberField(f);xx:=Evaluate(xp,a);yy:=Evaluate(yp,a);
 E:=FromXY(xx,yy);assert Discriminant(E) ne 0;P:=E![0,0];
 assert N*P eq E!0 and &and[(N div q)*P ne E!0:q in PrimeDivisors(N)];
 Q<U>:=PolynomialRing(K);bb2,bb4,bb6,bb8:=Explode(bInvariants(E));
 rt:=Roots(4*U^3+bb2*U^2+2*bb4*U+bb6);n2:=#rt;assert n2 in [0,1,3];
 if IsEven(N) then
  assert n2 gt 0;lower:=n2 eq 3 select 2*N else N;
  inv:=n2 eq 3 select [2,N] else [N];
 else
  lower:=n2 eq 3 select 4*N else (n2 eq 1 select 2*N else N);
  inv:=n2 eq 3 select [2,2*N] else (n2 eq 1 select [2*N] else [N]);
 end if;
 Note9(Sprintf("NONIC_EXACT_POINT index=%o N=%o field_degree=9 j_degree=%o rational_two_roots=%o torsion_lower=%o anticipated_invariants=%o time=%o",idx,N,Degree(MinimalPolynomial(jInvariant(E))),n2,lower,inv,Cputime(t0)));
 bound:=0;used:=0;df:=Discriminant(f);
 for p in PrimesUpTo(97) do
  if p eq 2 then continue;end if;
  if Numerator(df) mod p eq 0 or Denominator(df) mod p eq 0 then continue;end if;
  if &or[Denominator(z) mod p eq 0:z in Coefficients(f) cat Coefficients(xp) cat Coefficients(yp)] then continue;end if;
  fp:=Polynomial(GF(p),f);assert IsSquarefree(fp);
  for fg in Factorization(fp) do
   gg:=fg[1];d:=Degree(gg);
   if p^d gt 10^12 then continue;end if;
   if d eq 1 then L:=GF(p);aa:=-Coefficient(gg,0)/Coefficient(gg,1);
   else L<aa>:=ext<GF(p)|gg>;end if;
   xe:=EvalMod(xp,aa);ye:=EvalMod(yp,aa);
   if xe*ye eq 0 or xe^2*ye-xe eq 0 then continue;end if;
   rr:=(xe^2*ye-xe*ye+ye-1)/(xe^2*ye-xe);ss:=(xe*ye-ye+1)/(xe*ye);
   bb:=rr*ss*(rr-1);cc:=ss*(rr-1);
   b2:=(1-cc)^2-4*bb;b4:=-(1-cc)*bb;b6:=bb^2;b8:=-bb^3;
   if -b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6 eq 0 then continue;end if;
   Ee:=EllipticCurve([1-cc,-bb,-bb,0,0]);
    if Discriminant(Ee) eq 0 then continue;end if;
    PP:=Ee![0,0];assert N*PP eq Ee!0;
    nn:=#Ee;assert nn mod lower eq 0;
    bound:=GCD(bound,nn);used+:=1;
    Note9(Sprintf("NONIC_GOOD_REDUCTION index=%o p=%o residue_degree=%o order=%o gcd=%o time=%o",idx,p,d,nn,bound,Cputime(t0)));
   if bound eq lower then break;end if;
  end for;
  if bound eq lower then break;end if;
 end for;
 Note9(Sprintf("NONIC_TORSION_ORDER_BOUND index=%o N=%o lower=%o upper=%o primes_used=%o",idx,N,lower,bound,used));
 assert bound eq lower;
 Note9(Sprintf("NONIC_FULL_TORSION_CERTIFIED index=%o invariants=%o exact=true time=%o",idx,inv,Cputime(t0)));
end for;
Note9(Sprintf("NONIC_VANHOEIJ_EXAMPLES_DONE count=%o complete=true time=%o",#NonicExamples,Cputime(t0)));
delete out9;quit;
