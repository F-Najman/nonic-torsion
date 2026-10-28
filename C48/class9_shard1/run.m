SetColumns(0);SetQuitOnError(true);SetMemoryLimit(6*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "C48/data/portable_inputs.m";
load "common/pair_af_place_linear.m";
Y:=MDX11(2,24,GF(5));FC:=FunctionField(Curve(Y));F0,am:=AlgorithmicFunctionField(FC);F:=RationalExtensionRepresentation(F0);assert Genus(F) eq 13;B:=BaseField(F);k:=ConstantField(F);PP:=PolynomialRing(k);
function GetFunction(data)return F![B!(PP!v[1])/B!(PP!v[2]):v in data];end function;
assert Coefficients(DefiningPolynomial(F)) eq [B!(PP!v[1])/B!(PP!v[2]):v in PortableFieldEquation];
coords:=[F!am(z):z in Y`_coordinates];assert coords eq [GetFunction(v):v in PortableCoordinates];
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#PortableDegrees] do
 OO:=PortableFlags[j] eq 1 select OF else OI;I:=ideal<OO|[GetFunction(v):v in PortableIdealBases[j]]>;P:=Place(I);assert Degree(P) eq PortableDegrees[j];Append(~fb,P);
end for;assert #fb eq 810;
Z:=FreeAbelianGroup(9);G,pG:=quo<Z|[Z!v:v in PortableRelations]>;fbclasses:=[pG(Z!v):v in PortableClasses];A:=Divisor(fb[801]);ac:=fbclasses[801];assert Degree(A) eq 1;
H:=sub<G|[pG(Z!v):v in PortableRationalCusps]>;assert Invariants(H) eq [2,4,4,120,240];Q,q:=quo<G|[pG(Z!v):v in PortableRationalCusps] cat [ac]>;assert #Q eq 4608;

load "C48/portable_forward.m";
load "C48/data/parameter9/records.m";
load "C48/data/character_filter9/remaining.m";
assert #RemainingIndices eq 696;assert #Records eq 13513;
plane:=[F!am(FC.i):i in [1..2]];assert plane eq [F.1,F!B.1];
function UV(bc)
 b,c:=Explode(bc);t:=2*b/(1-c)-1;qq:=((1-c)*t^2-2)/(2*t);ss:=t*(qq^2-1)/(t*(qq^2+1)+2*qq);dd:=(ss+2)/qq;return [(ss+dd)/2,(ss-dd)/2];
end function;assert UV(coords) eq plane;
phi,witness:=C48PortableDirectRR(F,fb,fbclasses,G,A,"C48/logs/class9_shard1/identities.log");
R:=PolynomialRing(GF(5));K<z>:=ext<GF(5)|R!FieldPolynomial>;
ids:=[RemainingIndices[j]:j in [1..#RemainingIndices]|(j-1) mod 2 eq 1];out:=Open("C48/data/class9_shard1/results.py","w");
fprintf out,"shard=1\nexpected_indices=%o\nquotient_relations=%o\nquotient_images_of_G=%o\ntests=[\n",ids,[Eltseq(r):r in Rows(RelationMatrix(Q))],[Eltseq(q(pG(Z.i))):i in [1..9]];Flush(out);passing:=[];done:=0;
for i in ids do
 raw:=Records[i];bc:=[Evaluate(R!raw[j],z):j in [2,3]];uv:=UV(bc);P9:=C48PairAFPlaceLinear(plane,uv);assert Degree(P9) eq 9;
 a:=phi(Divisor(P9)-9*A);assert Eltseq(a)[9] eq 0;cl:=q(a);pass:=cl eq Q!0;if pass then Append(~passing,i);end if;
 fprintf out,"[%o,%o,%o,%o],\n",i,Eltseq(a),Eltseq(cl),pass select 1 else 0;Flush(out);done+:=1;
 if done mod 25 eq 0 or pass then printf "NONIC48_CLASS9_PROGRESS shard=1 done=%o passing=%o time=%o\n",done,#passing,Cputime(tt);end if;
end for;
fprintf witness,"]\n";delete witness;fprintf out,"]\npassing_indices=%o\ncomplete=True\n",passing;delete out;
printf "NONIC48_CLASS9_DONE shard=1 tested=%o passing=%o complete=true time=%o\n",done,#passing,Cputime(tt);quit;
