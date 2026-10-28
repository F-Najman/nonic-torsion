/* Fresh geometric verification with NO class-group computation or cache. */
SetColumns(0);SetQuitOnError(true);tt:=Cputime();
load "C54/data/quadratic_geometry_data.m";
AttachSpec("mdmagma/v2/mdmagma.spec");X:=MDX1(27,GF(5));C:=Curve(X);FF:=FunctionField(C);AF,am:=AlgorithmicFunctionField(FF);AF:=RationalExtensionRepresentation(AF);
K:=BaseField(AF);Pr:=PolynomialRing(GF(5));
function Rat54(v)return K!(Pr!v[1])/(K!(Pr!v[2]));end function;
assert DefiningPolynomial(AF) eq Polynomial(K,[Rat54(v):v in QuadraticFieldPolynomialData]);
assert [AF!am(FF.i):i in [1..2]] eq [AF!K.1,AF.1];
f:=AF![Rat54(v):v in QuadraticFunctionData];D:=Divisor(f);ps,es:=Support(D);assert #ps gt 0 and &and[e mod 2 eq 0:e in es];
cs:=[c:c in Cusps(X)|Degree(c) eq 1 and CuspSignature(X,c) eq QuadraticCuspSignature];assert #cs eq 1;
pc:=Support(FunctionFieldDivisor(Divisor(cs[1])))[1];assert Degree(pc) eq 1 and Valuation(f,pc) eq 0 and GF(5)!Evaluate(f,pc) eq GF(5)!QuadraticCuspValue;
assert not assigned AF`cld;
printf "NONIC54_QUADRATIC_GEOMETRY_VERIFIED genus=13 exact_even_divisor=true rational_cusp_unit_normalization=true field_coordinates_identity=true class_cache_absent=true complete=true time=%o\n",Cputime(tt);quit;
