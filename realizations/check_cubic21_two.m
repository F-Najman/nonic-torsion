SetColumns(0);SetQuitOnError(true);
R<x>:=PolynomialRing(Rationals());f:=x^3-3*x^2+3;assert IsIrreducible(f);
K<a>:=NumberField(f);y:=-a^2+2*a+1;r:=(a^2*y-a*y+y-1)/(a^2*y-a);s:=(a*y-y+1)/(a*y);b:=r*s*(r-1);c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]);P:=E![0,0];assert 21*P eq E!0 and 7*P ne E!0 and 3*P ne E!0;
S<u>:=PolynomialRing(K);b2,b4,b6,b8:=Explode(bInvariants(E));f2:=4*u^3+b2*u^2+2*b4*u+b6;
assert IsIrreducible(f2) and not IsSquare(Discriminant(f2));
printf "NONIC_CUBIC21_TWO_POINT_CHECK j=%o field_degree=3 two_cubic_irreducible=true two_discriminant_square=false full_two_splitting_relative_degree=6 complete=true\n",jInvariant(E);
quit;
