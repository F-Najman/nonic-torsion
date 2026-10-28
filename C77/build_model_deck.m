SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
CQ:=SmallModularCurve(11);jraw:=jInvariant(CQ,11);
CC:=ChangeRing(CQ,GF(2));assert IsNonsingular(CC) and Genus(CC) eq 1;
RC:=CoordinateRing(Ambient(CC));jn:=RC!Numerator(jraw);jd:=RC!Denominator(jraw);
K,kmap:=AlgorithmicFunctionField(FunctionField(CC));jK:=kmap(FunctionField(CC)!(jn/jd));
// The classical X0(7) parameter has j=(v^2+v+1)^4/v^7 in characteristic2.
Fv<v>:=FunctionField(GF(2));j7:=(v^2+v+1)^4/v^7;j7dual:=(v^2+v+1)^4/v;
Phi7:=ClassicalModularPolynomial(7);
assert Evaluate(Phi7,[j7,j7dual]) eq 0;
R<x>:=PolynomialRing(K);f:=x^8+jK*x^7+x^4+1;
assert IsIrreducible(f);AF<u>:=FunctionField(f);assert Genus(AF) eq 7;
jaf:=AF!jK;
cs:=Poles(jaf);widths:=[-Valuation(jaf,c):c in cs];
assert #cs eq 4 and &and[Degree(c) eq 1:c in cs] and Seqset(widths) eq {1,7,11,77};
printf "NONIC77_X0_MODEL genus7=true modular_polynomial_verified=true cusp_widths=%o time=%o\n",widths,Cputime(tt);
load "C77/AS_component.m";
RF<Y>:=PolynomialRing(AF);f2:=Y^2+Y+a2;assert IsIrreducible(f2);
// Direct cyclic cubic: the fine X1(7) family over its genus-zero parameter.
FB<t>:=FunctionField(GF(2));b7:=t^3-t^2;c7:=t^2-t;
E7u:=EllipticCurve([1-c7,-b7,-b7,FB!0,FB!0]);P7u:=E7u![0,0];assert Order(P7u) eq 7;
ju:=jInvariant(E7u);vu:=(t^2+t)/(t^3+t+1);
assert vu^8+ju*vu^7+vu^4+1 eq 0;
assert Evaluate(Numerator(ju),1/(t+1))/Evaluate(Denominator(ju),1/(t+1)) eq ju;
// j on X011 is j0(x)+j1(x)y, with y quadratic over F2(x).
jparts:=Eltseq(jK);assert #jparts eq 2 and jparts[2] ne 0;
kpoly:=DefiningPolynomial(K);assert Degree(kpoly) eq 2;
FX<X>:=FunctionField(FB);rx:=hom<BaseField(K)->FX|X>;
j0:=rx(jparts[1]);j1:=rx(jparts[2]);aa:=rx(Coefficient(kpoly,1));bb:=rx(Coefficient(kpoly,0));
normrel:=((FX!ju)+j0)^2+aa*j1*((FX!ju)+j0)+bb*j1^2;
f3:=Numerator(normrel);assert Degree(f3) eq 12 and IsIrreducible(f3);
E3<s>:=FunctionField(f3);
kx:=hom<BaseField(K)->E3|s>;
yval:=((E3!ju)+kx(jparts[1]))/kx(jparts[2]);
kemb:=hom<K->E3|kx,yval>;
afemb:=hom<AF->E3|kemb,E3!vu>;
assert afemb(jaf) eq E3!ju;
tau:=hom<E3->E3|hom<FB->E3|1/((E3!t)+1)>,s>;
assert tau(tau(tau(E3!t))) eq E3!t and tau(tau(tau(s))) eq s;
assert &and[tau(afemb(v)) eq afemb(v):v in [AF.1,AF!K.1,AF!K!BaseField(K).1]];
assert Genus(E3) eq 19;
printf "NONIC77_DIRECT_CUBIC_VERIFIED degree12=true genus19=true exact_base_embedding=true deck_order3=true time=%o\n",Cputime(tt);
QY<z>:=PolynomialRing(E3);fy:=QY![afemb(c):c in Coefficients(f2)];assert IsIrreducible(fy);
CY<zz>:=FunctionField(fy);jY:=CY!E3!ju;
assert Genus(CY) eq 37;
printf "NONIC77_DIRECT_GENUS37_VERIFIED time=%o\n",Cputime(tt);
cys:=Poles(jY);cywidths:=[-Valuation(jY,P):P in cys];
assert {*<cywidths[i],Degree(cys[i])>:i in [1..#cys]*} eq {*<77,1>^^6,<7,2>^^3,<11,3>^^2,<1,6>^^1*};
assert &+[cywidths[i]*Degree(cys[i]):i in [1..#cys]] eq 576;
printf "NONIC77_GLOBAL_CUSPS_VERIFIED signature=%o width_degree_pairs=%o jdegree576=true time=%o\n",[Degree(P):P in cys],[<cywidths[i],Degree(cys[i])>:i in [1..#cys]],Cputime(tt);
print "NONIC77_DIRECT_GLOBAL_GENUS37_MODEL_DONE";

 SetColumns(0);dt:=Cputime();assert Genus(CY) eq 37;
 assert &+[cywidths[i]*Degree(cys[i]):i in [1..#cys]] eq 576;
 cs77:=[cys[i]:i in [1..#cys]|cywidths[i] eq 77];assert #cs77 eq 6;
 assert BaseField(CY) cmpeq E3;
 deck2:=hom<CY->CY|map<E3->CY|a:->CY!a>,CY.1+1>;
 deck3:=hom<CY->CY|map<E3->CY|a:->CY!tau(a)>,CY.1>;
 generatorsY:=[CY.1,CY!E3.1,CY!FB.1];
 for a in generatorsY do
  assert deck2(deck2(a)) eq a and deck3(deck3(deck3(a))) eq a;
  assert deck2(deck3(a)) eq deck3(deck2(a));
 end for;
 labels77:=[];
 for P in cs77 do
  tv:=Valuation(CY!FB.1,P) lt 0 select 2 else Integers()!Evaluate(CY!FB.1,P);
  assert tv in [0,1,2] and Valuation(CY.1,P) ge 0;
  zv:=Integers()!Evaluate(CY.1,P);assert zv in [0,1];
  Append(~labels77,<tv,zv>);
 end for;
 assert Seqset(labels77) eq {<i,j>:i in [0,1,2],j in [0,1]};
 pp2:=[Index(labels77,<q[1],1-q[2]>):q in labels77];
 pp3:=[Index(labels77,<q[1] eq 0 select 1 else (q[1] eq 1 select 2 else 0),q[2]>):q in labels77];
 printf "NONIC77_DECK_LABELS labels=%o perm2=%o perm3=%o time=%o\n",labels77,pp2,pp3,Cputime(dt);
 Omega77,om77:=SpaceOfHolomorphicDifferentials(CY);B77:=[om77(v):v in Basis(Omega77)];
 CM77:=Matrix(GF(2),[Eltseq(Cartier(w)@@om77):w in B77]);
 U77:=RowSpace(CM77*(CM77-IdentityMatrix(GF(2),37)));assert Dimension(U77) eq 32;
 eta77:=Differential(CY!FB.1);hs77:=[w/eta77:w in B77];MM77:=[];
 for k in [1,2] do
  sigma:=k eq 1 select deck2 else deck3;de:=Differential(sigma(CY!FB.1));rows:=[];
  for j in [1..37] do
   w:=sigma(hs77[j])*de;Append(~rows,Eltseq(w@@om77));
   if j mod 5 eq 0 then printf "NONIC77_DECK_DIFFERENTIAL_PROGRESS deck=%o basis=%o/37 time=%o\n",k,j,Cputime(dt);end if;
  end for;
  Append(~MM77,Matrix(GF(2),rows));
 end for;
 M2:=MM77[1];M3:=MM77[2];I37:=IdentityMatrix(GF(2),37);
 assert M2^2 eq I37 and M3^3 eq I37 and M2*M3 eq M3*M2;
 assert Rank(M2) eq 37 and Rank(M3) eq 37;
 assert M2*CM77 eq CM77*M2 and M3*CM77 eq CM77*M3;
 assert #Seqset([M2^a*M3^b:a in [0,1],b in [0,1,2]]) eq 6;
 assert &and[v*M in U77:v in Basis(U77),M in [M2,M3]];
 printf "NONIC77_DECK_DIFFERENTIALS_VERIFIED group_order6=true commutes_Cartier=true image_dimension32=true time=%o\n",Cputime(dt);

print "NONIC77_MODEL_AND_DECK_READY";

load "C77/rr9_continue.m";
