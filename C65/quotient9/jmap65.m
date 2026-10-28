SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
load "C65/data/quotient9/canonical_data.m";
gg:=9;prec:=650;pp:=2;F5:=GF(2);
P:=ProjectiveSpace(F5,gg-1);R:=CoordinateRing(P);mons:=[R.(z[1])*R.(z[2]):z in pairs];eqs:=[&+[v[i]*mons[i]:i in [1..#pairs]]:v in Rows(QuadRelations)];C:=Curve(P,eqs);IC:=Ideal(C);GB:=GroebnerBasis(IC);LM:=[LeadingMonomial(f):f in GB];
function StdMons(dd)
 ss:={R!1};for k in [1..dd] do cand:={m*R.i:m in ss,i in [1..gg]};ss:={m:m in cand|not exists{l:l in LM|IsDivisibleBy(m,l)}};assert #ss eq (k eq 1 select gg else (2*k-1)*(gg-1));end for;return Sort(Setseq(ss));
end function;
PS<q>:=PowerSeriesRing(F5,prec);Q:=[PS!Eltseq(v)+O(q^prec):v in Rows(Submatrix(qcoeffs,1,1,gg,prec))];
function ExpMat(ms)return Matrix(F5,[[Coefficient(f,n):n in [0..prec-1]] where f:=Evaluate(m,Q):m in ms]);end function;
Delta:=q+O(q^prec);
for n in [1..prec-1] do Delta*:= (1-q^n)^24;end for;
found:=false;
for aa in [3..8] do
 Sa:=StdMons(aa);Sm:=StdMons(aa+6);Ea:=ExpMat(Sa);Em:=ExpMat(Sm);assert Rank(Ea) eq #Sa and Rank(Em) eq #Sm;KK:=Kernel(VerticalJoin(Ea,Em));
 printf "NONIC65_JMAP_E4_SECTION a=%o dimension=%o time=%o\n",aa,Dimension(KK),Cputime(tt);
 if Dimension(KK) eq 0 then continue;end if;
 V:=BasisMatrix(KK);VG:=Submatrix(V,1,1,Nrows(V),#Sa);DG:=Matrix(F5,[[Coefficient(Delta*(PS!Eltseq(v)+O(q^prec)),n):n in [0..prec-1]]:v in Rows(VG*Ea)]);
 WW:=Kernel(VerticalJoin(DG,Em));
 printf "NONIC65_JMAP_BOTH_SECTION a=%o dimension=%o time=%o\n",aa,Dimension(WW),Cputime(tt);
 if Dimension(WW) gt 0 then
  rr:=WW.1;cv:=Vector(F5,Eltseq(rr)[1..Nrows(V)]);vv:=cv*V;gv:=Vector(F5,Eltseq(vv)[1..#Sa]);av:=-Vector(F5,Eltseq(vv)[#Sa+1..#Sa+#Sm]);bv:=-Vector(F5,Eltseq(rr)[Nrows(V)+1..Nrows(V)+#Sm]);
  assert gv ne 0 and gv*Ea eq av*Em;
  G:=PS!Eltseq(gv*Ea)+O(q^prec);assert bv*Em eq Vector(F5,[Coefficient(Delta*G,n):n in [0..prec-1]]);
  found:=true;aDegree:=aa;break;
 end if;
end for;assert found;
assert prec gt 48*(aDegree+6);
PrintFile("C65/data/quotient9/jmap_q_data.m",Sprintf("N:=65;genus:=9;p:=2;precision:=%o;SturmBound:=%o;SectionDegree:=%o;Gdegree:=%o;StdGExponents:=%o;StdJExponents:=%o;Gvec:=%o;NumeratorVec:=%o;DenominatorVec:=%o;",prec,48*(aDegree+6),aDegree+6,aDegree,[Exponents(m):m in Sa],[Exponents(m):m in Sm],Eltseq(gv),Eltseq(av),Eltseq(bv)):Overwrite:=true);
printf "NONIC65_JMAP_Q_CERTIFIED E4cubed_mod2_equals1=true weight=%o Sturm=%o precision=%o time=%o\n",2*(aDegree+6),48*(aDegree+6),prec,Cputime(tt);
load "C65/data/quotient9/plane_model_data.m";
K<t>:=FunctionField(F5);PF:=PolynomialRing(F5);RY<y>:=PolynomialRing(K);
function DecK(c)return (K!PF!c[1])/(K!PF!c[2]);end function;
f:=RY![DecK(c):c in AFPolynomialData];AF<g>:=FunctionField(f);assert Degree(AF) eq AF_degree and Genus(AF) eq 9;
function DecAF(cs)return &+[AF!DecK(cs[i])*g^(i-1):i in [1..#cs]];end function;
if AFBaseCoordinate eq 1 then uAF:=AF!t;vAF:=g;else uAF:=g;vAF:=AF!t;end if;
Z:=[DecAF(cs):cs in CanonicalCoordinateData];
fv:=AF!0;for i in [1..#PlaneExponents] do e:=PlaneExponents[i];c:=PlaneCoefficients[i];if e[2] gt 0 then fv+:=(AF!F5!(c*e[2]))*uAF^e[1]*vAF^(e[2]-1);end if;end for;
assert fv ne 0;Z:=[fv*w:w in Z];DD:=LCM([Denominator(c):c in Eltseq(w),w in Z]);Z:=[AF!DD*w:w in Z];assert &and[Denominator(c) eq 1:c in Eltseq(w),w in Z];
pows:=AssociativeArray();for i in [1..gg] do pows[<i,0>]:=AF!1;for k in [1..aDegree+6] do pows[<i,k>]:=pows[<i,k-1>]*Z[i];end for;end for;
ms:=[&*[pows[<i,e[i]>]:i in [1..gg]]:e in [Exponents(m):m in Sm]];
j:=(&+[av[i]*ms[i]:i in [1..#Sm]])/(&+[bv[i]*ms[i]:i in [1..#Sm]]);assert Degree(j) eq 168;
cs:=Poles(j);cuspdata:=Sort([<Degree(c),-Valuation(j,c)>:c in cs]);assert &+[Degree(c)*(-Valuation(j,c)):c in cs] eq 168;

function CodeK(w)return <[Integers()!a:a in Eltseq(Numerator(w))],[Integers()!a:a in Eltseq(Denominator(w))]>;end function;
function CodeAF(w)return [CodeK(c):c in Eltseq(w)];end function;
PrintFile("C65/data/quotient9/finite_j_model_data.m",Sprintf("N:=65;genus:=9;p:=2;AF_degree:=%o;AFBaseCoordinate:=%o;AFPolynomialData:=%o;JFunctionData:=%o;CanonicalCoordinateData:=%o;CuspDegreeWidth:=%o;",Degree(AF),AFBaseCoordinate,AFPolynomialData,CodeAF(j),[CodeAF(w):w in Z],cuspdata):Overwrite:=true);
printf "NONIC65_FINITE_J_MODEL_READY jdegree168=true cusp_degree_width=%o time=%o\n",cuspdata,Cputime(tt);
HH:=[u:u in[1..64]|GCD(u,65)eq 1 and KroneckerSymbol(65,u)eq 1];S:=CuspidalSubspace(ModularSymbolsH(65,HH,2,1));I:=IdentityMatrix(Rationals(),9);ords:=[];
for q in [2,3,7,11,17,19,23,29,31,37,41,43] do
 T:=HeckeOperator(S,q);D:=(T^2-HeckeOperator(S,q^2))/q;v:=Integers()!Determinant(I-T+q*D);Append(~ords,<q,v>);
 printf "NONIC65_RATIONAL_ORDER_BOUND prime=%o order=%o gcd_so_far=%o time=%o\n",q,v,GCD([o[2]:o in ords]),Cputime(tt);
end for;
PrintFile("C65/data/quotient9/rational_order_bounds.m",Sprintf("N:=65;genus:=9;Orders:=%o;RationalTorsionOrderBound:=%o;",ords,GCD([o[2]:o in ords])):Overwrite:=true);
print "NONIC65_JMODEL_AND_BOUNDS_COMPLETE";quit;
