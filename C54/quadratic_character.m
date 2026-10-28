/* Exact unramified quadratic evaluation obstruction, after PROVEN Pic(F5).
   Inputs AF5,G5,phi,coords,X5,C5; no class queries are used for evaluations. */
assert #ConstantField(AF5) eq 5 and coords eq [AF5!BaseField(AF5).1,AF5.1];
GC,qcinv,qcphi:=ClassGroup(AF5:Proof:=true);assert GC cmpeq G5;
ii:=Invariants(GC);assert ii eq [3,3,342175446,0];aa2:=(ii[3] div 2)*GC.3;assert Order(aa2) eq 2;
function ReduceZero54(D)
 E,r,A:=Reduction(D);answer:=E+r*A;assert Degree(answer) eq 0;return answer;
end function;
baseDiv:=ReduceZero54(qcinv(GC.3));DD2:=DivisorGroup(AF5)!0;
for bit in IntegerToSequence(ii[3] div 2,2) do
 if bit eq 1 then DD2:=ReduceZero54(DD2+baseDiv);end if;
 baseDiv:=ReduceZero54(2*baseDiv);
end for;
assert Degree(DD2) eq 0 and qcphi(DD2) eq aa2;
bb:=Basis(2*DD2:Reduction:=false,Simplification:="None");assert #bb eq 1;quadf:=bb[1];
OF:=MaximalOrderFinite(AF5);OI:=MaximalOrderInfinite(AF5);ID,JD:=Ideals(2*DD2);
assert ID*(quadf*OF) eq ideal<OF|1> and JD*(quadf*OI) eq ideal<OI|1>;
DF:=Divisor(quadf);ps,es:=Support(DF);assert #ps gt 0 and &and[e mod 2 eq 0:e in es];
found:=false;
for cc in Cusps(X5) do
 if Degree(cc) ne 1 then continue;end if;
 DC:=FunctionFieldDivisor(Divisor(cc));pcc:=Support(DC)[1];
 if Valuation(quadf,pcc) eq 0 then basecusp:=cc;fc:=GF(5)!Evaluate(quadf,pcc);found:=true;break;end if;
end for;
assert found and fc ne 0;quadcoeff:=Eltseq(quadf);
function QuadraticCharacterXY(xy)
 K:=Parent(xy[1]);xx,yy:=Explode(xy);value:=K!0;power:=K!1;
 for c in quadcoeff do
  de:=Evaluate(Denominator(c),xx);if de eq 0 then return 0;end if;
  value+:=(Evaluate(Numerator(c),xx)/de)*power;power*:=yy;
 end for;
 if value eq 0 then return 0;end if;
 sign:=(value/(K!fc))^((#K-1) div 2);assert sign in {K!1,-K!1};return sign eq 1 select 1 else -1;
end function;
function EncRat54(a)return [[Integers()!c:c in Coefficients(Numerator(a))],[Integers()!c:c in Coefficients(Denominator(a))]];end function;
if not assigned QuadraticFile then QuadraticFile:="C54/data/quadratic_character_data.json";end if;
outq:=Open(QuadraticFile,"w");
fprintf outq,"{\"field_polynomial\":%o,\"filter_function\":%o,\"base_cusp_signature\":%o,\"base_cusp_value\":%o,\"function_divisor_even\":true,\"both_orders_exact\":true,\"Pic_invariants\":%o,\"H5_odd_order\":472473,\"complete\":true}\n",[EncRat54(c):c in Coefficients(DefiningPolynomial(AF5))],[EncRat54(c):c in quadcoeff],CuspSignature(X5,basecusp),Integers()!fc,ii;delete outq;
printf "NONIC54_QUADRATIC_CHARACTER_READY exact_even_divisor=true both_orders_exact=true base_cusp_signature=%o base_value=%o function_degree=%o time=%o\n",CuspSignature(X5,basecusp),fc,Degree(Numerator(DF)),Cputime(tt);
