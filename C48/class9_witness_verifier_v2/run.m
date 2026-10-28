SetColumns(0);SetQuitOnError(true);SetMemoryLimit(5*10^9);tt:=Cputime();AttachSpec("mdmagma/v2/mdmagma.spec");load "C48/data/portable_inputs.m";load "C48/data/parameter9/records.m";load "C48/data/class9_witness_verifier_v2/inputs.m";
Y:=MDX11(2,24,GF(5));FC:=FunctionField(Curve(Y));F0,am:=AlgorithmicFunctionField(FC);F:=RationalExtensionRepresentation(F0);B:=BaseField(F);kk:=ConstantField(F);R:=PolynomialRing(kk);assert Genus(F) eq 13;
function Rat(v)return B!(R!v[1])/B!(R!v[2]);end function;assert Coefficients(DefiningPolynomial(F)) eq [Rat(v):v in PortableFieldEquation];
function GetFunction(data)return F![Rat(v):v in data];end function;
plane:=[F!am(FC.i):i in [1..2]];assert plane eq [F.1,F!B.1];coords:=[F!am(v):v in Y`_coordinates];assert coords eq [GetFunction(v):v in PortableCoordinates];
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[* *];
for j in [1..810] do OO:=PortableFlags[j] eq 1 select OF else OI;I:=ideal<OO|[GetFunction(v):v in PortableIdealBases[j]]>;assert Degree(Place(I)) eq PortableDegrees[j];Append(~fb,I);end for;
function DivisorIdeals(sparse)
 IF:=ideal<OF|1>;II:=ideal<OI|1>;degree:=0;seen:={Integers()|};cls:=[0:i in [1..9]];
 for row in sparse do j,e:=Explode(row);assert j in [1..810] and e gt 0 and j notin seen;Include(~seen,j);degree+:=e*PortableDegrees[j];cls:=[cls[k]+e*PortableClasses[j][k]:k in [1..9]];
  if PortableFlags[j] eq 1 then IF*:=fb[j]^e;else II*:=fb[j]^e;end if;
 end for;return IF,II,degree,cls;
end function;
AF,AI,ad,ac:=DivisorIdeals([[801,1]]);assert ad eq 1;
function PairDescriptor(values)
 K:=Parent(values[1]);kk:=PrimeField(K);ee:=LCM([Degree(MinimalPolynomial(v,kk)):v in values]);found:=false;
 for bits in [[1,0],[0,1]] cat [[1,a]:a in [1..4]] do gen:=&+[K!bits[j]*values[j]:j in [1,2]];pol:=MinimalPolynomial(gen,kk);if Degree(pol) eq ee then chosen:=bits;found:=true;break;end if;end for;assert found;
 M:=Matrix(kk,[Eltseq(gen^j,kk):j in [0..ee-1]]);assert Rank(M) eq ee;cs:=[Eltseq(Solution(M,Vector(kk,Eltseq(v,kk)))):v in values];assert &and[&+[K!cs[j][i+1]*gen^i:i in [0..ee-1]] eq values[j]:j in [1,2]];
 return [ee] cat chosen cat [Integers()!u:u in Coefficients(pol)] cat [Integers()!u:u in &cat cs];
end function;
function UV(b,c)
 t:=2*b/(1-c)-1;q:=((1-c)*t^2-2)/(2*t);ss:=t*(q^2-1)/(t*(q^2+1)+2*q);dd:=(ss+2)/q;return [(ss+dd)/2,(ss-dd)/2];
end function;
K<z>:=ext<GF(5)|R!FieldPolynomial>;mods:=[PortableRelations[j][j]:j in [1..8]];done:=0;
for rec in ReplayCases do
 idx,cls,w:=Explode(rec);assert w[2] eq 0 and w[3] eq cls and [[Integers()!v:v in row]:row in w[6]] eq [[801,1]];EF,EI,ed,ec:=DivisorIdeals(w[4]);RF,RI,rd,rc:=DivisorIdeals(w[5]);assert ed eq 13 and rd eq 13;
 aa:=GetFunction(w[7]);assert aa ne 0;ID:=EF/(RF*(aa*OF));JD:=EI/(RI*(aa*OI));IP:=ID*AF^9;JP:=JD*AI^9;assert IsIntegral(IP) and IsIntegral(JP);
 factors:=[*<v[1],v[2]>:v in Factorization(IP)*];for v in Factorization(JP) do Append(~factors,<v[1],v[2]>);end for;assert #factors eq 1 and factors[1][2] eq 1;pl:=Place(factors[1][1]);assert Degree(pl) eq 9;
 PF,PI:=Ideals(Divisor(pl));assert PF eq IP and PI eq JP;DF,DI:=Ideals(Divisor(pl)-9*Divisor(Place(PortableFlags[801] eq 1 select AF else AI)));assert DF eq ID and DI eq JD;
 actual:=[Evaluate(f,pl):f in plane];raw:=Records[idx];b:=Evaluate(R!raw[2],z);c:=Evaluate(R!raw[3],z);expected:=UV(b,c);assert PairDescriptor(actual) eq PairDescriptor(expected);
 assert ec[9]-rc[9] eq cls[9] and &and[(ec[j]-rc[j]-cls[j]) mod mods[j] eq 0:j in [1..8]];done+:=1;
 if done mod 50 eq 0 then printf "NONIC48_CLASS9_WITNESS_PROGRESS checked=%o time=%o\n",done,Cputime(tt);end if;
end for;
assert done eq 696 and not assigned F`cld;
printf "NONIC48_CLASS9_WITNESS_VERIFIER_DONE checked696=true both_orders_exact=true actual_point_bindings=true class_coordinates_recomputed=true no_RR_or_class_search=true complete=true time=%o\n",Cputime(tt);quit;
