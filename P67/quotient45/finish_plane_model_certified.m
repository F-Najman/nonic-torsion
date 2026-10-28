// Prepared only: no native launch without root resource allocation.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
load "P67/data/quotient45/canonical_data.m";
load "P67/data/quotient45/jmap_q_data.m";
load "P67/data/quotient45/plane_differential_certified_data.m";
load "P67/data/quotient45/canonical_change_basis_certified.m";
K<t>:=FunctionField(GF(2));PF:=PolynomialRing(GF(2));RY<y>:=PolynomialRing(K);
function DecK(c)return (K!PF!c[1])/(K!PF!c[2]);end function;
f:=RY![DecK(c):c in AFPolynomialData];AF<g>:=FunctionField(f);assert Degree(AF) eq AF_degree and Genus(AF) eq 45;
function DecAF(cs)return &+[AF!DecK(cs[i])*g^(i-1):i in [1..#cs]];end function;
if AFBaseCoordinate eq 1 then uAF:=AF!t;vAF:=g;else uAF:=g;vAF:=AF!t;end if;
assert &+[uAF^e[1]*vAF^e[2]:e in PlaneExponents] eq 0;
BC:=[DecAF(cs):cs in DifferentialCoefficientData];assert #BC eq 45;
// Plane adjoints have the common differential denominator F_v.
// Remove this common projective factor before high-degree substitutions.
fv:=&+[uAF^e[1]*vAF^(e[2]-1):e in PlaneExponents|IsOdd(e[2])];assert fv ne 0;
BC:=[fv*w:w in BC];Z:=[&+[CanonicalChangeBasis[i,j]*BC[j]:j in [1..45]]:i in [1..45]];
denoms:=[Denominator(c):c in Eltseq(w),w in Z];DD:=LCM(denoms);Z:=[AF!DD*w:w in Z];
assert &and[Denominator(c) eq 1:c in Eltseq(w),w in Z];
printf "NONIC67H3_ADJOINT_COORDINATES_READY genus=45 AF_degree=%o maximum_coefficient_degree=%o time=%o\n",Degree(AF),Max([Degree(Numerator(c)):c in Eltseq(w),w in Z]),Cputime(tt);
cz:=[Z[p[1]]*Z[p[2]]:p in pairs];assert &and[&+[v[i]*cz[i]:i in [1..1035]] eq 0:v in Rows(QuadRelations)];
la0:=&+[ProjectionRows[1,i]*Z[i]:i in [1..45]];
assert (&+[ProjectionRows[2,i]*Z[i]:i in [1..45]])/la0 eq uAF and (&+[ProjectionRows[3,i]*Z[i]:i in [1..45]])/la0 eq vAF;
printf "NONIC67H3_CANONICAL_PULLBACK_EXACT all903quadrics=true projection_inverse=true time=%o\n",Cputime(tt);
pows:=AssociativeArray();for i in [1..45] do pows[<i,0>]:=AF!1;for k in [1..10] do pows[<i,k>]:=pows[<i,k-1>]*Z[i];end for;end for;
ms8:=[&*[pows[<i,e[i]>]:i in [1..45]]:e in Std10Exponents];
num:=&+[NumeratorVec[i]*ms8[i]:i in [1..836]];den:=&+[DenominatorVec[i]*ms8[i]:i in [1..836]];j:=num/den;
assert Degree(j) eq 748;cs:=Poles(j);assert #cs eq 12 and Sort([Degree(pl):pl in cs]) eq ([1:i in[1..11]]cat[11]) and Sort([-Valuation(j,pl):pl in cs]) eq ([1:i in[1..11]]cat[67]);
printf "NONIC67H3_PLANE_JMAP_EXACT genus=45 jdegree=748 rational_cusps=11 nonrational_cusp_degree=11 pole_orders=1^11_67 time=%o\n",Cputime(tt);
SZ:=[&+[Diamond2Matrix[i,j]*Z[j]:j in [1..45]]:i in [1..45]];la:=&+[ProjectionRows[1,i]*SZ[i]:i in [1..45]];su:=(&+[ProjectionRows[2,i]*SZ[i]:i in [1..45]])/la;sv:=(&+[ProjectionRows[3,i]*SZ[i]:i in [1..45]])/la;
// The canonical diamond matrix comes from exact q-eigenform characters.
// It preserves the modular j-function; finite place keys use that matrix directly.
function CodeK(w)return <[Integers()!a:a in Eltseq(Numerator(w))],[Integers()!a:a in Eltseq(Denominator(w))]>;end function;
function CodeAF(w)return [CodeK(c):c in Eltseq(w)];end function;
out:="P67/data/quotient45/plane_model_certified_data.m";
PrintFile(out,Sprintf("N:=67;genus:=45;plane_degree:=%o;AF_degree:=%o;AFBaseCoordinate:=%o;AFGeneratorCoordinate:=%o;PlaneExponents:=%o;ProjectionRows:=Matrix(GF(2),3,45,%o);AFPolynomialData:=%o;CanonicalCoordinateData:=%o;JFunctionData:=%o;SigmaUData:=%o;SigmaVData:=%o;Diamond2Matrix:=Matrix(GF(2),45,45,%o);",plane_degree,Degree(AF),AFBaseCoordinate,AFGeneratorCoordinate,PlaneExponents,Eltseq(ProjectionRows),[CodeK(w):w in Eltseq(DefiningPolynomial(AF))],[CodeAF(w):w in Z],CodeAF(j),CodeAF(su),CodeAF(sv),Eltseq(Diamond2Matrix)):Overwrite:=true);
PrintFile(out,"print \"NONIC67H3_LOW_DEGREE_MODEL_DATA_READY\";");printf "NONIC67H3_LOW_DEGREE_MODEL_CERTIFIED genus=45 AF_degree=%o j_degree=748 diamond2_modular_matrix_order=11 time=%o\n",Degree(AF),Cputime(tt);quit;
