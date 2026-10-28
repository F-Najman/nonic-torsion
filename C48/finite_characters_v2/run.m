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
S:=sub<Q|>;selected:=[];
for j in [1..810] do
 if PortableDegrees[j] lt 3 then continue;end if;
 z:=q(fbclasses[j]);if z notin S then Append(~selected,j);S:=sub<Q|[Q!v:v in Generators(S)] cat [z]>;end if;
 if #S eq #Q then break;end if;
end for;assert #S eq #Q;
ZZ:=FreeAbelianGroup(#selected);ff:=hom<ZZ->Q|[q(fbclasses[j]):j in selected]>;assert IsSurjective(ff);
preimages:=[Eltseq(Q.i@@ff):i in [1..Ngens(Q)]];
function Bucket(P)return [Degree(P),IsFinite(P) select 1 else 0] cat (IsFinite(P) select [Integers()!v:v in Coefficients(Minimum(P))] else []);end function;
buckets:=AssociativeArray();for j in [1..810] do ky:=Bucket(fb[j]);if IsDefined(buckets,ky) then Append(~buckets[ky],j);else buckets[ky]:=[j];end if;end for;
function FBIndex(P)
 ky:=Bucket(P);assert IsDefined(buckets,ky);hits:=[j:j in buckets[ky]|fb[j] eq P];assert #hits eq 1;return hits[1];
end function;

plane:=[F!am(FC.i):i in [1..2]];
function UVFromBC(bc)
 b,c:=Explode(bc);t:=2*b/(1-c)-1;q:=((1-c)*t^2-2)/(2*t);ss:=t*(q^2-1)/(t*(q^2+1)+2*q);dd:=(ss+2)/q;return [(ss+dd)/2,(ss-dd)/2];
end function;
assert UVFromBC(coords) eq plane;print "NONIC48_EXPLICIT_TARGET_INVERSE_VERIFIED";

labels:=[<eps,eta,a>:eps,eta in [0,1],a in [1,5,7,11]];actions:=[* *];imageindices:=[];
for label in labels do
 eps,eta,a:=Explode(label);ims:=[];inds:=[];
 for j in selected do
  b,c:=Explode([Evaluate(v,fb[j]):v in coords]);K:=Parent(b);E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P2:=E![0,0];Q24:=E![b,b];assert 24*Q24 eq E!0 and 12*Q24 ne E!0 and 8*Q24 ne E!0;
  bc:=DerickxNormalForm_bc(E,P2+eps*(12*Q24),a*Q24+eta*P2);pl:=C48PairAFPlaceLinear(plane,UVFromBC(bc));assert Degree(pl) eq PortableDegrees[j];idx:=FBIndex(pl);Append(~inds,idx);Append(~ims,q(fbclasses[idx]));
 end for;
 action:=hom<Q->Q|[&+[preimages[i][j]*ims[j]:j in [1..#selected]]:i in [1..Ngens(Q)]]>;assert #Kernel(action) eq 1 and IsSurjective(action);Append(~actions,action);Append(~imageindices,inds);
 print "NONIC48_AUT_ACTION",label,Cputime(tt);
end for;
out:=Open("C48/data/finite_characters_v2/action_data.py","w");
fprintf out,"quotient_relations=%o\nquotient_images_of_G_generators=%o\nselected_factorbase_indices=%o\nselected_preimages=%o\nlabels=%o\nimages_of_selected=%o\naction_images=%o\ncomplete=True\n",[Eltseq(r):r in Rows(RelationMatrix(Q))],[Eltseq(q(pG(Z.i))):i in [1..9]],selected,preimages,[[l[i]:i in [1..3]]:l in labels],imageindices,[[Eltseq(action(Q.i)):i in [1..Ngens(Q)]]:action in actions];delete out;

// All eight finite-Jacobian2-torsion functions, using portable geometry only.
T,qT:=quo<G|[ac]>;assert #T eq 4246732800;SS:=sub<T|>;chosen:=[];
for j in [1..810] do
 if PortableDegrees[j] lt 3 then continue;end if;
 z:=qT(fbclasses[j]);if z notin SS then Append(~chosen,j);SS:=sub<T|[T!v:v in Generators(SS)] cat [z]>;end if;if #SS eq #T then break;end if;
end for;assert #SS eq #T;
ST:=FreeAbelianGroup(#chosen);ft:=hom<ST->T|[qT(fbclasses[j]):j in chosen]>;assert IsSurjective(ft);
function RZero(D)
 E,r,Ared:=Reduction(D);ans:=E+r*Ared;assert Degree(ans) eq 0;return ans;
end function;
VS:=VectorSpace(GF(2),8);span:=sub<VS|>;basisindices:=[];basisvectors:=[];
for j in [1..810] do
 if PortableDegrees[j] lt 3 then continue;end if;vv:=Eltseq(fbclasses[j]-PortableDegrees[j]*ac);assert vv[9] eq 0;v:=VS!vv[1..8];if v notin span then Append(~basisindices,j);Append(~basisvectors,v);span:=sub<VS|basisvectors>;end if;if Dimension(span) eq 8 then break;end if;
end for;assert #basisindices eq 8;MM:=Matrix(basisvectors);basepi:=LocalUniformizer(fb[801]);assert Valuation(basepi,fb[801]) eq 1;
function RatData(z)return [[Integers()!v:v in Coefficients(Numerator(z))],[Integers()!v:v in Coefficients(Denominator(z))]];end function;
functions:=[* *];charrows:=[];basevals:=[];inverseVectors:=[];
for i in [1..Ngens(T)] do
 ni:=Order(T.i);assert ni mod 2 eq 0;v:=Eltseq(((ni div 2)*T.i)@@ft);Append(~inverseVectors,v);DD:=DivisorGroup(F)!0;
 printf "NONIC48_FINITE2_FUNCTION_START i=%o small_generators=%o maxcoeff=%o time=%o\n",i,#chosen,Max([Abs(z):z in v]),Cputime(tt);
 for j in [1..#chosen] do
  n:=v[j];if n eq 0 then continue;end if;bb:=Divisor(fb[chosen[j]])-PortableDegrees[chosen[j]]*A;if n lt 0 then bb:=-bb;n:=-n;end if;
  while n gt 0 do if n mod 2 eq 1 then DD:=RZero(DD+bb);end if;n div:=2;if n gt 0 then bb:=RZero(2*bb);end if;end while;
 end for;
 assert Degree(DD) eq 0;bs:=Basis(2*DD:Reduction:=false,Simplification:="None");assert #bs eq 1;f:=bs[1];ID,JD:=Ideals(2*DD);assert ID*(f*OF) eq ideal<OF|1> and JD*(f*OI) eq ideal<OI|1>;
 Df:=Divisor(f);ps,es:=Support(Df);assert &and[e mod 2 eq 0:e in es];vf:=Valuation(f,fb[801]);fc:=GF(5)!Evaluate(f/basepi^vf,fb[801]);assert fc ne 0;
 signs:=[];for j in basisindices do
  vv:=Valuation(f,fb[j]);assert vv mod 2 eq 0;
  if vv eq 0 then value:=Evaluate(f,fb[j]);else pj:=LocalUniformizer(fb[j]);value:=Evaluate(f/pj^vv,fb[j]);end if;
  K:=Parent(value);sg:=(value/K!fc)^((#K-1) div 2);assert sg in {K!1,-K!1};Append(~signs,sg eq 1 select 0 else 1);
 end for;
 cr:=Solution(Transpose(MM),Vector(GF(2),signs));Append(~charrows,cr);Append(~basevals,Integers()!fc);Append(~functions,[RatData(z):z in Eltseq(f)]);
 printf "NONIC48_FINITE2_FUNCTION_DONE i=%o char=%o time=%o\n",i,Eltseq(cr),Cputime(tt);
end for;
CM:=Matrix(charrows);assert Rank(CM) eq 8;HM:=Matrix(GF(2),[[Integers()!x:x in v[1..8]]:v in PortableRationalCusps]);kill:=Kernel(CM*Transpose(HM));KK:=Matrix([v*CM:v in Basis(kill)]);assert Rank(KK) eq 4;
out:=Open("C48/data/finite_characters_v2/function_data.py","w");fprintf out,"functions=%o\nbase_values=%o\ncharacter_rows=%o\nannihilating_combinations=%o\nbasis_indices=%o\nbasis_vectors=%o\nchosen_factorbase_indices=%o\ninverse_vectors=%o\ncomplete=True\n",functions,basevals,[[Integers()!v:v in Eltseq(r)]:r in charrows],[[Integers()!v:v in Eltseq(r)]:r in Basis(kill)],basisindices,[[Integers()!v:v in Eltseq(r)]:r in basisvectors],chosen,inverseVectors;delete out;
printf "NONIC48_FINITE_CHARACTERS_DONE actions16=true full_character_rank=8 H_annihilating_rank=4 no_classgroup_search=true complete=true time=%o\n",Cputime(tt);quit;
