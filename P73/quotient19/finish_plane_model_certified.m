SetColumns(0);SetQuitOnError(true);SetMemoryLimit(48*10^9);tt:=Cputime();
load "P73/data/quotient19/canonical_data.m";
load "P73/data/quotient19/jmap_q_data.m";
load "P73/data/quotient19/plane_differential_certified_data.m";
load "P73/data/quotient19/canonical_change_basis_certified.m";
K<t>:=FunctionField(GF(2));PF:=PolynomialRing(GF(2));RY<y>:=PolynomialRing(K);
function DecK(c)return (K!PF!c[1])/(K!PF!c[2]);end function;
f:=RY![DecK(c):c in AFPolynomialData];AF<g>:=FunctionField(f);assert Degree(AF) eq AF_degree and Genus(AF) eq 19;
function DecAF(cs)return &+[AF!DecK(cs[i])*g^(i-1):i in [1..#cs]];end function;
if AFBaseCoordinate eq 1 then uAF:=AF!t;vAF:=g;else uAF:=g;vAF:=AF!t;end if;
assert &+[uAF^e[1]*vAF^e[2]:e in PlaneExponents] eq 0;
BC:=[DecAF(cs):cs in DifferentialCoefficientData];assert #BC eq 19;
// Plane adjoints have the common differential denominator F_v.
// Remove this common projective factor before high-degree substitutions.
fv:=&+[uAF^e[1]*vAF^(e[2]-1):e in PlaneExponents|IsOdd(e[2])];assert fv ne 0;
BC:=[fv*w:w in BC];Z:=[&+[CanonicalChangeBasis[i,j]*BC[j]:j in [1..19]]:i in [1..19]];
denoms:=[Denominator(c):c in Eltseq(w),w in Z];DD:=LCM(denoms);Z:=[AF!DD*w:w in Z];
assert &and[Denominator(c) eq 1:c in Eltseq(w),w in Z];
printf "NONIC73_ADJOINT_COORDINATES_READY genus=19 AF_degree=10 maximum_coefficient_degree=%o time=%o\n",Max([Degree(Numerator(c)):c in Eltseq(w),w in Z]),Cputime(tt);
cz:=[Z[p[1]]*Z[p[2]]:p in pairs];assert &and[&+[v[i]*cz[i]:i in [1..190]] eq 0:v in Rows(QuadRelations)];
la0:=&+[ProjectionRows[1,i]*Z[i]:i in [1..19]];
assert (&+[ProjectionRows[2,i]*Z[i]:i in [1..19]])/la0 eq uAF and (&+[ProjectionRows[3,i]*Z[i]:i in [1..19]])/la0 eq vAF;
printf "NONIC73_CANONICAL_PULLBACK_EXACT all136quadrics=true projection_inverse=true time=%o\n",Cputime(tt);
pows:=AssociativeArray();for i in [1..19] do pows[<i,0>]:=AF!1;for k in [1..9] do pows[<i,k>]:=pows[<i,k-1>]*Z[i];end for;end for;
ms8:=[&*[pows[<i,e[i]>]:i in [1..19]]:e in Std9Exponents];
num:=&+[NumeratorVec[i]*ms8[i]:i in [1..306]];den:=&+[DenominatorVec[i]*ms8[i]:i in [1..306]];j:=num/den;
assert Degree(j) eq 296;cs:=Poles(j);assert #cs eq 8 and &and[Degree(pl) eq 1:pl in cs]and Sort([-Valuation(j,pl):pl in cs]) eq [1,1,1,1,73,73,73,73];
printf "NONIC73_PLANE_JMAP_EXACT genus=19 jdegree=296 rational_cusps=8 pole_orders=1^4_73^4 time=%o\n",Cputime(tt);
SZ:=[&+[Diamond5Matrix[i,j]*Z[j]:j in [1..19]]:i in [1..19]];la:=&+[ProjectionRows[1,i]*SZ[i]:i in [1..19]];su:=(&+[ProjectionRows[2,i]*SZ[i]:i in [1..19]])/la;sv:=(&+[ProjectionRows[3,i]*SZ[i]:i in [1..19]])/la;
// The canonical diamond matrix comes from exact q-eigenform characters.
// It preserves the modular j-function; finite place keys use that matrix directly.
function CodeK(w)return <[Integers()!a:a in Eltseq(Numerator(w))],[Integers()!a:a in Eltseq(Denominator(w))]>;end function;
function CodeAF(w)return [CodeK(c):c in Eltseq(w)];end function;
out:="P73/data/quotient19/plane_model_certified_data.m";
PrintFile(out,Sprintf("N:=73;genus:=19;plane_degree:=%o;AF_degree:=%o;AFBaseCoordinate:=%o;AFGeneratorCoordinate:=%o;PlaneExponents:=%o;ProjectionRows:=Matrix(GF(2),3,19,%o);AFPolynomialData:=%o;CanonicalCoordinateData:=%o;JFunctionData:=%o;SigmaUData:=%o;SigmaVData:=%o;Diamond5Matrix:=Matrix(GF(2),19,19,%o);",plane_degree,Degree(AF),AFBaseCoordinate,AFGeneratorCoordinate,PlaneExponents,Eltseq(ProjectionRows),[CodeK(w):w in Eltseq(DefiningPolynomial(AF))],[CodeAF(w):w in Z],CodeAF(j),CodeAF(su),CodeAF(sv),Eltseq(Diamond5Matrix)):Overwrite:=true);
PrintFile(out,"print \"NONIC73_LOW_DEGREE_MODEL_DATA_READY\";");printf "NONIC73_LOW_DEGREE_MODEL_CERTIFIED genus=19 AF_degree=10 j_degree=296 diamond5_modular_matrix_order=4 time=%o\n",Cputime(tt);quit;
