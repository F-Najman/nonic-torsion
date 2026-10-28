SetColumns(0);SetQuitOnError(true);SetMemoryLimit(6*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "C48/degree8/two_isogeny_48_to_2_24.m";
F<t>:=FunctionField(Rationals());aa:=2*t;bb:=t^2-1;cc:=t^2+1;
E:=EllipticCurve([F!0,aa^4+bb^4,F!0,aa^4*bb^4,F!0]);Q8:=E![aa*bb*(aa+cc)*(bb+cc),aa*bb*cc*(aa+bb)*(aa+cc)*(bb+cc)];P2:=E![-aa^4,0];assert 8*Q8 eq E!0 and 4*Q8 ne E!0 and P2 ne 4*Q8;
jj:=jInvariant(E);assert Max(Degree(Numerator(jj)),Degree(Denominator(jj))) eq 48;print "NONIC48_GENUS0_GENERIC_VERIFIED degree_j=48 order8=true independent2=true";
function TOrbit(t)
 K:=Parent(t);if t in {K!0,K!1,-K!1} then return [t];end if;
 u:=(t-1)/(t+1);v:=(t+1)/(t-1);return [t,-t,1/t,-1/t,u,-u,v,-v];
end function;
function Key(bc,d)
 return Minimum([Eltseq(bc[1]^(5^i),GF(5)) cat Eltseq(bc[2]^(5^i),GF(5)):i in [0..d-1]]);
end function;
for dd in [3,4,5] do
 K<z>:=GF(5^dd);RR<x>:=PolynomialRing(K);seen:={K|};targetKeys:={};sourceKeys:={};autoKeys:={};fieldpoly:=Eltseq(DefiningPolynomial(K,GF(5)));targetOrbits:=0;curveCases:=0;noncuspt:=0;
 for t in K do
  if t in seen then continue;end if;
  torb:=TOrbit(t);for u in torb,i in [0..dd-1] do Include(~seen,u^(5^i));end for;
  aa:=2*t;bb:=t^2-1;cc:=t^2+1;A:=aa^4+bb^4;B:=aa^4*bb^4;if B*(A^2-4*B) eq 0 then continue;end if;
  E:=EllipticCurve([K!0,A,K!0,B,K!0]);Q8:=E![aa*bb*(aa+cc)*(bb+cc),aa*bb*cc*(aa+bb)*(aa+cc)*(bb+cc)];P2:=E![-aa^4,0];assert 8*Q8 eq E!0 and 4*Q8 ne E!0 and 2*P2 eq E!0 and P2 ne E!0 and P2 ne 4*Q8;curveCases+:=1;
  for rt in Roots(DivisionPolynomial(E,3)) do
   xr:=rt[1];yes,yr:=IsSquare(xr^3+A*xr^2+B*xr);if not yes then continue;end if;assert yr ne 0;R3:=E![xr,yr];assert R3 ne E!0 and 3*R3 eq E!0;Q24:=3*Q8+R3;assert 24*Q24 eq E!0 and 12*Q24 ne E!0 and 8*Q24 ne E!0;
   pairs:=[];keys:=[];
   for eps,eta in [0,1],a in [1,5,7,11] do
    bc:=DerickxNormalForm_bc(E,P2+eps*(12*Q24),a*Q24+eta*P2);Append(~pairs,bc);Append(~keys,Key(bc,dd));
   end for;
   deg:=LCM([Degree(MinimalPolynomial(v,GF(5))):v in pairs[1]]);if deg ne dd then continue;end if;
   canonical:=Minimum(keys);if canonical in autoKeys then continue;end if;Include(~autoKeys,canonical);targetOrbits+:=1;
   localkeys:={};
   for k in [1..16] do
    ky:=keys[k];if ky in localkeys then continue;end if;Include(~localkeys,ky);assert ky notin targetKeys;Include(~targetKeys,ky);b,c:=Explode(pairs[k]);yes,sb:=IsSquare(b);if not yes then continue;end if;
    d:=b*(1-b-c);Ap:=-2*c;Bp:=c^2-4*d;ES:=EllipticCurve([K!0,Ap,K!0,Bp,K!0]);
    for sg in [K!1,-K!1] do
     xx:=c+2*b+2*sg*sb;yy:=8*b*xx^2/(Bp-xx^2);PS:=ES![xx,yy];assert 48*PS eq ES!0 and 24*PS ne ES!0 and 16*PS ne ES!0;
     EB,LB:=C48TwoIsogenousModuli(ES,PS);back:=DerickxNormalForm_bc(EB,LB`P,LB`Q);assert Key(back,dd) eq ky;
     sk:=Key(TateNormalForm_bc(ES,PS),dd);assert sk notin sourceKeys;Include(~sourceKeys,sk);
    end for;
   end for;
  end for;
 end for;
 assert #seen eq #K;
 expectedTarget:=[16,96,672][dd-2];expectedSource:=[16,64,656][dd-2];assert #targetKeys eq expectedTarget and #sourceKeys eq expectedSource;
 printf "NONIC48_GENUS0_FINITE_CONTROL degree=%o curves=%o target_Aut16_orbits=%o targetclosed=%o sourceclosed=%o all_complete=true time=%o\n",dd,curveCases,targetOrbits,#targetKeys,#sourceKeys,Cputime(tt);
end for;
printf "NONIC48_PARAMETER_PILOT_DONE generic_j_degree48=true all_small_source_and_target_Newton_counts_matched=true complete=true time=%o\n",Cputime(tt);quit;
