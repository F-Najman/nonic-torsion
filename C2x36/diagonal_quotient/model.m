/* Experimental quotient by diamond13, not an exclusion certificate. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Y:=MDX11(2,12,GF(5));C:=Curve(Y);FF:=FunctionField(C);AF,am:=AlgorithmicFunctionField(FF);assert Genus(AF) eq 1;
E:=EllipticCurve([am(FF!a):a in Y`_E]);P:=E![am(FF!a):a in Y`_P];Q:=E![am(FF!a):a in Y`_Q];assert 12*Q eq E!0 and 6*Q ne E!0 and 4*Q ne E!0 and P ne 6*Q;
T:=4*Q;R<x>:=PolynomialRing(AF);EP,ph:=IsogenyFromKernel(E,x-T[1]);psi:=DivisionPolynomial(EP,3);fac:=Factorization(psi);print "C236_DELTA13_QUOTIENT_CUBIC_FACTORS",[Degree(h[1]):h in fac],Cputime(tt);assert Sort([Degree(h[1]):h in fac]) eq [1,3];cub:=[h[1]:h in fac|Degree(h[1]) eq 3][1];
Z<root>:=FunctionField(cub);assert Degree(Z) eq 3;assert Genus(Z) eq 9;zcoords:=[Z!am(FF!f):f in Y`_coordinates] cat [root];zj:=Z!jInvariant(E);
