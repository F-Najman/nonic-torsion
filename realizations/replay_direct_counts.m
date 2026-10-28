/* Independent finite-curve order verification by summing quadratic characters.
   This uses no elliptic-curve point-count algorithm or cached finite orders. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t0:=Cputime();
load "common/data/inputs/vanhoeij_nonic_input.m";
load "realizations/data/selected_reductions.m";
R<T>:=PolynomialRing(Rationals());checked:=0;
function Ev9(f,a)
 return &+[Parent(a)!Numerator(Coefficient(f,i))/(Parent(a)!Denominator(Coefficient(f,i)))*a^i:i in [0..Degree(f)]];
end function;
for entry in SelectedReductions do
 idx,rr:=Explode(entry);N,fc,xc,yc:=Explode(NonicExamples[idx]);
 f:=R!fc;xp:=R!xc;yp:=R!yc;bound:=0;
 for rec in rr do
  p,d,expected:=Explode(rec);assert IsPrime(p) and p ne 2;
  assert Numerator(Discriminant(f)) mod p ne 0 and Denominator(Discriminant(f)) mod p ne 0;
  assert &and[Denominator(z) mod p ne 0:z in Coefficients(f) cat Coefficients(xp) cat Coefficients(yp)];
  fg:=[g[1]:g in Factorization(Polynomial(GF(p),f))|Degree(g[1]) eq d];assert #fg gt 0;
  // Check every place of the selected residue degree: the producer may have
  // used any one of them. At least one must have the recorded order.
  orders:=[];
  for g in fg do
   if d eq 1 then L:=GF(p);aa:=-Coefficient(g,0)/Coefficient(g,1);
   else L<aa>:=ext<GF(p)|g>;end if;
   x0:=Ev9(xp,aa);y0:=Ev9(yp,aa);
   if x0*y0 eq 0 or x0^2*y0-x0 eq 0 then continue;end if;
   r:=(x0^2*y0-x0*y0+y0-1)/(x0^2*y0-x0);s:=(x0*y0-y0+1)/(x0*y0);
   b:=r*s*(r-1);c:=s*(r-1);
   b2:=(1-c)^2-4*b;b4:=-(1-c)*b;b6:=b^2;b8:=-b^3;
   if -b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6 eq 0 then continue;end if;
   E:=EllipticCurve([1-c,-b,-b,0,0]);
   if Discriminant(E) eq 0 then continue;end if;
   a1,a2,a3,a4,a6:=Explode(aInvariants(E));count:=1;
   for z in L do
    disc:=(a1*z+a3)^2+4*(z^3+a2*z^2+a4*z+a6);
    if disc eq 0 then count+:=1;
    elif disc^((#L-1) div 2) eq 1 then count+:=2;
    end if;
   end for;
   Append(~orders,count);
  end for;
  assert expected in orders;bound:=GCD(bound,expected);checked+:=1;
  printf "NONIC_DIRECT_COUNT index=%o p=%o residue_degree=%o expected=%o observed=%o exact=true\n",idx,p,d,expected,orders;
 end for;
 assert bound eq N;
end for;
printf "NONIC_DIRECT_COUNTS_COMPLETED examples=19 reductions=%o time=%o complete=true\n",checked,Cputime(t0);quit;
