/* Direct C4 cover of the genus11 quotient; new degree9 section test. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");load "common/pair_af_place_linear.m";
first:=eval Read("C2x36/data/filter_quotient18_linear_arithmetic.txt");row:=[r:r in first[3]|r[1] eq 8 and r[2] eq 1601][1];def:=[r[2]:r in first[2]|r[1] eq 8][1];R5:=PolynomialRing(GF(5));L<zz>:=ext<GF(5)|R5!def>;key:=row[3];bc:=[Evaluate(R5!key[1..8],zz),Evaluate(R5!key[9..16],zz)];EL:=EllipticCurve([L|0,bc[2],0,bc[1]*(1-bc[1]-bc[2]),0]);PL:=EL![0,0];QL:=EL![bc[1],bc[1]];
A<t>:=FunctionField(GF(5));r:=t^2-t+1;c9:=t*(r-1);b9:=r*c9;E9:=EllipticCurve([A|1-c9,-b9,-b9,0,0]);j9:=jInvariant(E9);RW<W>:=PolynomialRing(A);polyW:=16*(W^2+16)^3-j9*W^4;assert IsIrreducible(polyW);Z<w>:=FunctionField(polyW);assert Genus(Z) eq 11;
B9,C9:=Explode(TateNormalForm_bc(EL,4*QL));tv:=C9^2/(B9-C9);wv:=4*((9*QL)[1]-(18*QL)[1])/(PL[1]-(PL+18*QL)[1]);q2:=C48PairAFPlaceLinear([Z!t,w],[tv,wv]);assert Degree(q2) eq 2 and Valuation(Z!j9,q2) ge 0;assert Evaluate(Numerator(j9),tv)/Evaluate(Denominator(j9),tv) eq jInvariant(EL);
print "C236_D9_QUARTIC_BASE_POINT_READY",Cputime(tt);
aa:=aInvariants(E9);a2:=Z!(aa[2]+aa[1]^2/4);a4:=Z!(aa[4]+aa[1]*aa[3]/2);a6:=Z!(aa[5]+aa[3]^2/4);
RZ<X>:=PolynomialRing(Z);cubic:=X^3+a2*X^2+a4*X+a6;vA:=a2+3*X;vB:=a4+2*a2*X+3*X^2;rel:=w^2*vA^2-4*(w^2+4)*vB;linear:=GCD(cubic,rel);assert Degree(linear) eq 1;rho:=-Coefficient(linear,0)/Coefficient(linear,1);
alpha:=a2+3*rho;beta:=a4+2*a2*rho+3*rho^2;quartic:=X^4-2*alpha*beta*X^2+alpha^2*beta^2-4*beta^3;assert IsIrreducible(quartic);S<y4>:=FunctionField(quartic);assert Degree(S) eq 4;
x4shift:=(y4^2-S!(alpha*beta))/(2*S!beta);assert x4shift^2 eq S!beta;delta:=4*x4shift/(S!w);assert delta^2 eq S!(alpha^2-4*beta);
ES:=EllipticCurve([S|0,a2,0,a4,a6]);PS:=ES![S!rho+(-S!alpha+delta)/2,0];P9:=ES![0,S!(aa[3]/2)];S4:=ES![S!rho+x4shift,y4];assert 2*S4 eq ES![S!rho,0] and 4*S4 eq ES!0 and 9*P9 eq ES!0 and 3*P9 ne ES!0;QS:=S4+7*P9;assert 36*QS eq ES!0 and PS ne 18*QS;
print "C236_D9_QUARTIC_SOURCE_CURVE_READY",Cputime(tt);

/* Character decomposition avoids normalizing the genus41 source. */
iota:=GF(5)!2;assert iota^2 eq -1;kap:=4*beta/w;zeta:=y4+(S!iota)*(S!kap)*x4shift/y4;fk:=16*beta^3*(1+2*iota/w)^2;assert zeta^4 eq S!fk;

 assert Valuation(fk,q2) eq 0;K2,mp2:=ResidueClassField(q2);vals:=[Evaluate(z,q2):z in [Z!t,w,rho,fk]];poly2:=[Integers()!c:c in Coefficients(DefiningPolynomial(K2,GF(5)))];coeffs:=[[Integers()!c:c in Eltseq(z,GF(5))]:z in vals];Write("C2x36/data/orbit_section_criterion/inert_basepoint_values.txt",Sprint([*poly2,coeffs*]):Overwrite:=true);print "C236_D9_EXPLICIT_INERT_BASEPOINT_VALUES_PASS complete=true",Cputime(tt);quit;
