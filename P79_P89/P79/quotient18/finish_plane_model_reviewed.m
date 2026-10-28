SetColumns(0);SetQuitOnError(true);SetMemoryLimit(48*10^9);tt:=Cputime();
load "P79_P89/data/P79/quotient18/canonical_data.m";
load "P79_P89/data/P79/quotient18/jmap_q_data.m";
load "P79_P89/data/P79/quotient18/plane_differential_certified_data.m";
load "P79_P89/data/P79/quotient18/canonical_change_basis_certified.m";
K<t>:=FunctionField(GF(2));PF:=PolynomialRing(GF(2));RY<y>:=PolynomialRing(K);
function DecK(c)return (K!PF!c[1])/(K!PF!c[2]);end function;
f:=RY![DecK(c):c in AFPolynomialData];AF<g>:=FunctionField(f);assert Degree(AF) eq AF_degree and Genus(AF) eq 18;
function DecAF(cs)return &+[AF!DecK(cs[i])*g^(i-1):i in [1..#cs]];end function;
if AFBaseCoordinate eq 1 then uAF:=AF!t;vAF:=g;else uAF:=g;vAF:=AF!t;end if;
assert &+[uAF^e[1]*vAF^e[2]:e in PlaneExponents] eq 0;
BC:=[DecAF(cs):cs in DifferentialCoefficientData];assert #BC eq 18;
// Plane adjoints have the common differential denominator F_v.
// Remove this common projective factor before high-degree substitutions.
fv:=&+[uAF^e[1]*vAF^(e[2]-1):e in PlaneExponents|IsOdd(e[2])];assert fv ne 0;
BC:=[fv*w:w in BC];Z:=[&+[CanonicalChangeBasis[i,j]*BC[j]:j in [1..18]]:i in [1..18]];
denoms:=[Denominator(c):c in Eltseq(w),w in Z];DD:=LCM(denoms);Z:=[AF!DD*w:w in Z];
assert &and[Denominator(c) eq 1:c in Eltseq(w),w in Z];
printf "NONIC79_ADJOINT_COORDINATES_READY genus=18 AF_degree=13 maximum_coefficient_degree=%o time=%o\n",Max([Degree(Numerator(c)):c in Eltseq(w),w in Z]),Cputime(tt);
cz:=[Z[p[1]]*Z[p[2]]:p in pairs];assert &and[&+[v[i]*cz[i]:i in [1..171]] eq 0:v in Rows(QuadRelations)];
la0:=&+[ProjectionRows[1,i]*Z[i]:i in [1..18]];
assert (&+[ProjectionRows[2,i]*Z[i]:i in [1..18]])/la0 eq uAF and (&+[ProjectionRows[3,i]*Z[i]:i in [1..18]])/la0 eq vAF;
printf "NONIC79_CANONICAL_PULLBACK_EXACT all120quadrics=true projection_inverse=true time=%o\n",Cputime(tt);
pows:=AssociativeArray();for i in [1..18] do pows[<i,0>]:=AF!1;for k in [1..8] do pows[<i,k>]:=pows[<i,k-1>]*Z[i];end for;end for;
ms8:=[&*[pows[<i,e[i]>]:i in [1..18]]:e in Std8Exponents];
num:=&+[NumeratorVec[i]*ms8[i]:i in [1..255]];den:=&+[DenominatorVec[i]*ms8[i]:i in [1..255]];j:=num/den;
assert Degree(j) eq 240;cs:=Poles(j);assert #cs eq 4 and Sort([Degree(pl):pl in cs]) eq [1,1,1,3] and Sort([-Valuation(j,pl):pl in cs]) eq [1,1,1,79];
printf "NONIC79_PLANE_JMAP_EXACT genus=18 jdegree=240 cusp_degrees=1^3_3 pole_orders=1^3_79 time=%o\n",Cputime(tt);
SZ:=[&+[Diamond3Matrix[i,j]*Z[j]:j in [1..18]]:i in [1..18]];la:=&+[ProjectionRows[1,i]*SZ[i]:i in [1..18]];su:=(&+[ProjectionRows[2,i]*SZ[i]:i in [1..18]])/la;sv:=(&+[ProjectionRows[3,i]*SZ[i]:i in [1..18]])/la;
// The canonical diamond matrix comes from exact q-eigenform characters.
// It preserves the modular j-function; finite place keys use that matrix directly.
function CodeK(w)return <[Integers()!a:a in Eltseq(Numerator(w))],[Integers()!a:a in Eltseq(Denominator(w))]>;end function;
function CodeAF(w)return [CodeK(c):c in Eltseq(w)];end function;
out:="P79_P89/data/P79/quotient18/plane_model_reviewed_data.m";
PrintFile(out,Sprintf("N:=79;genus:=18;plane_degree:=%o;AF_degree:=%o;AFBaseCoordinate:=%o;AFGeneratorCoordinate:=%o;PlaneExponents:=%o;ProjectionRows:=Matrix(GF(2),3,18,%o);AFPolynomialData:=%o;CanonicalCoordinateData:=%o;JFunctionData:=%o;SigmaUData:=%o;SigmaVData:=%o;Diamond3Matrix:=Matrix(GF(2),18,18,%o);",plane_degree,Degree(AF),AFBaseCoordinate,AFGeneratorCoordinate,PlaneExponents,Eltseq(ProjectionRows),[CodeK(w):w in Eltseq(DefiningPolynomial(AF))],[CodeAF(w):w in Z],CodeAF(j),CodeAF(su),CodeAF(sv),Eltseq(Diamond3Matrix)):Overwrite:=true);
PrintFile(out,"print \"NONIC79_LOW_DEGREE_MODEL_DATA_READY\";");printf "NONIC79_LOW_DEGREE_MODEL_CERTIFIED genus=18 AF_degree=13 j_degree=240 diamond3_modular_matrix_order=3 time=%o\n",Cputime(tt);quit;
