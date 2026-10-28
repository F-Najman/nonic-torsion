/* Model of Y = X_1(75)/<26> over F_2 as a degree-20 cover of P^1. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
R<T,x,y>:=PolynomialRing(GF(2),3);B:=FieldOfFractions(R);
F25:=eval Read("mdmagma/models_X1_n/FFFc25.txt");
xx:=B!x;yy:=B!y;r:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));s:=(xx*yy-yy+1)/(xx*yy);
b:=r*s*(r-1);c:=s*(r-1);
// The exact degree5 modular unit from the archived vanHoeij gonality table.
f3:=b;f4:=c;f5:=b-c;f6:=s-1;f7:=s-r;f8:=r*s-2*r+1;f9:=s^2-s-r+1;
f10:=xx-yy+1;f12:=xx-yy;f13:=xx^3*yy-xx^2*yy^2-xx^2*yy+xx*yy^2-yy+1;
h:=f3*f4*f6^3*f7^2*f10*f12^3/(f5*f8^2*f9^2*f13^3);
N:=R!Numerator(h);D:=R!Denominator(h);
print "N75_GON5_CANONICAL_RESULTANT_START",TotalDegree(N),TotalDegree(D),Cputime(tt);
res:=Resultant(F25,N-T*D,3);assert res ne 0;
K<t>:=RationalFunctionField(GF(2));S<X>:=PolynomialRing(K);proj:=hom<R->S|[S!t,X,S!0]>;
pol:=proj(res);factors:=Factorization(pol);
print "N75_GON5_CANONICAL_RESULTANT_FACTORS",[<Degree(z[1]),z[2]>:z in factors],Cputime(tt);
fs:=[z[1]:z in factors|Degree(z[1]) eq 5];assert #fs eq 1;
L<a>:=FunctionField(fs[1]);SY<Y>:=PolynomialRing(L);ev:=hom<R->SY|[SY!t,SY!a,Y]>;
gg:=GCD(ev(F25),ev(N-T*D));assert Degree(gg) eq 1;
ya:=-Coefficient(gg,0)/Coefficient(gg,1);assert Evaluate(ev(F25),ya) eq 0;
xa:=a;ra:=(xa^2*ya-xa*ya+ya-1)/(xa*(xa*ya-1));sa:=(xa*ya-ya+1)/(xa*ya);
ba:=ra*sa*(ra-1);ca:=sa*(ra-1);E:=EllipticCurve([1-ca,-ba,-ba,L!0,L!0]);P:=E![0,0];
assert 25*P eq E!0 and 5*P ne E!0 and Genus(L) eq 12;
print "N75_GON5_CANONICAL_BASE_MODEL_VERIFIED",Degree(L),Genus(L),Cputime(tt);
S3<W3>:=PolynomialRing(L);j25:=jInvariant(E);f:=W3^4+W3^3+1/j25;assert Degree(f) eq 4;
assert IsIrreducible(f);M<w>:=FunctionField(f);
print "N75_GON5_CANONICAL_DEGREE20_GENUS_START",Cputime(tt);
assert Genus(M) eq 73;
print "N75_GON5_CANONICAL_GENUS73_MODEL_VERIFIED",Cputime(tt);

function Gon5Value(xx,yy)
 r:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));s:=(xx*yy-yy+1)/(xx*yy);
 b:=r*s*(r-1);c:=s*(r-1);
 return b*c*(s-1)^3*(s-r)^2*(xx-yy+1)*(xx-yy)^3/((b-c)*(r*s-2*r+1)^2*(s^2-s-r+1)^2*(xx^3*yy-xx^2*yy^2-xx^2*yy+xx*yy^2-yy+1)^3);
end function;
AttachSpec("mdmagma/v2/mdmagma.spec");
xy2:=TateNormalForm_xy(E,2*P);h2:=Gon5Value(xy2[1],xy2[2]);
sig0:=hom<K->L|h2>;sig:=hom<L->L|sig0,xy2[1]>;
assert sig(ya) eq xy2[2] and sig(j25) eq j25;
assert &and[(sig^10)(x) eq x:x in [L!t,a]];
dg:=hom<M->M|map<L->M|x:->M!sig(x)>,w>;
assert &and[(dg^10)(x) eq x:x in [M!t,M!a,w]];
print "N75_D2_AUTOMORPHISM_VERIFIED",Cputime(tt);
bcs:=Poles(j25);assert Sort([Degree(p):p in bcs]) eq [1:i in [1..10]] cat [4,4,10];
selected:=[[i:i in [1..#bcs]|Degree(bcs[i]) eq d][1]:d in [1,4]];
sourcecs:=[];
for i in selected do
 p:=bcs[i];ff,gg:=TwoGenerators(p);u:=UniformizingElement(p);
 raw:=Decomposition(M,p);above:=[q:q in raw|Valuation(M!ff,q) gt 0 and Valuation(M!gg,q) gt 0];
 assert #above eq 2 and Sort([Valuation(M!u,q):q in above]) eq [1,3];
 assert &and[Degree(q) eq Degree(p):q in above];sourcecs cat:=above;
 print "N75_ORBIT_REPRESENTATIVE_CUSPS",i,Degree(p),[<-Valuation(M!j25,q),Degree(q)>:q in above],Cputime(tt);
end for;
assert #sourcecs eq 4;
Om,om:=SpaceOfHolomorphicDifferentials(M);bs:=[om(v):v in Basis(Om)];
C:=Matrix(GF(2),[Eltseq(Cartier(omega)@@om):omega in bs]);
I:=IdentityMatrix(GF(2),73);U:=RowSpace(C-I);BM:=BasisMatrix(U);

assert Genus(M) eq 73 and Dimension(U) eq 71;
q75:=[q:q in sourcecs|Degree(q) eq 1 and -Valuation(M!j25,q) eq 75][1];
V,vm:=DifferentialSpace(9*Divisor(q75));W:=sub<VectorSpace(GF(2),73)|[Eltseq(vm(v)@@om):v in Basis(V)]>;
r9:=Dimension(U)-Dimension(U meet W);assert r9 eq 9;
printf "NONIC75_PURE9_RR_COMPLETE genus73=true image71=true width75=true rank=%o covered10=true time=%o\n",r9,Cputime(tt);quit;
