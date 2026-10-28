SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
load "C33/data/F5_rr/input33.m";
k:=GF(5);K<t>:=FunctionField(k);PK:=PolynomialRing(k);RY<y>:=PolynomialRing(K);
function Rat(v)return (K!PK!v[1])/(K!PK!v[2]);end function;
F<u>:=FunctionField(RY![Rat(v):v in FieldEquation]);assert Degree(F)eq 12 and Genus(F)eq 21;
function Fun(v)return F![Rat(w):w in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);blocks:=[* *];
for z in GlobalCuspBlocks do
 sig,deg,pp:=Explode(z);DD:=DivisorGroup(F)!0;
 for v in pp do dd,m,flag,ib:=Explode(v);O:=flag eq 1 select OF else OI;II:=ideal<O|[Fun(b):b in ib]>;assert IsPrime(II);P:=Place(II);assert Degree(P)eq dd;DD+:=m*Divisor(P);end for;
 assert Degree(DD)eq deg;Append(~blocks,DD);
end for;

load "C33/data/F5_rr/basis33.m";
load "C33/data/F5_rr/jet_tables33.m";
A:=&+[ACoefficients[j]*blocks[j]:j in[1..17]];assert Degree(A)eq 48;fs:=[Fun(v):v in BasisFunctions];assert #fs eq 28;
for f in fs do assert f ne 0;for P in Poles(f)do assert Valuation(f,P)+Valuation(A,P)ge 0;end for;end for;
dn:=LCM([LCM([Denominator(c):c in Eltseq(f)]):f in fs]);md:=Max([Max([Degree(PK!(dn*c)):c in Eltseq(f)]):f in fs]);
IM:=Matrix(k,[&cat[[Coefficient(PK!(dn*c),n):n in[0..md]]:c in Eltseq(f)]:f in fs]);assert Rank(IM)eq 28;
PrintFile("C33/data/F5_rr/coefficient_independence33.json",Sprintf("{\"rows\":28,\"columns\":%o,\"max_degree\":%o,\"common_denominator\":%o,\"rank\":28}",Ncols(IM),md,[Integers()!c:c in Coefficients(dn)]):Overwrite:=true);
checks:=0;
for rec in JetTables do
 j,k0,dd,ord,orders,M:=Explode(rec);ps,es:=Support(blocks[j]);P:=ps[k0];assert Degree(P)eq dd and ord eq Valuation(A,P);
 maxb:=Max(orders);LS,comp:=Completion(F,P:Precision:=100);z:=LS.1;vv:=[z^ord*comp(f):f in fs];assert &and[Valuation(v)ge 0 and AbsolutePrecision(v)ge maxb:v in vv];
 N:=Matrix(k,[&cat[[Integers()!c:c in Eltseq(Coefficient(v,n),k)]:n in[0..maxb-1]]:v in vv]);
 for m in Seqset(orders)do
  if m eq 0 then continue;end if;SM:=Submatrix(M,1,1,28,m*dd);NM:=Submatrix(N,1,1,28,m*dd);W:=BasisMatrix(Kernel(SM));assert W*NM eq 0 and Rank(SM)eq Rank(NM);checks+:=1;
 end for;
 printf "NONIC33_F5_FRESH_LOCAL_KERNEL_VERIFIED block=%o place=%o degree=%o orders=%o time=%o\n",j,k0,dd,Sort(Setseq(Seqset(orders))),Cputime(tt);
end for;
for i in[1..10]do
 J:=ZeroMatrix(k,28,0);for rec in JetTables do if rec[5][i]gt 0 then J:=HorizontalJoin(J,Submatrix(rec[6],1,1,28,rec[3]*rec[5][i]));end if;end for;assert Ncols(J)eq 39 and Rank(J)eq 28;
end for;assert not assigned F`cld;
printf "NONIC33_F5_FRESH_POLE_JET_VERIFIER_COMPLETE basis28_poles_valid=true coefficient_independence28=true all_prefix_kernels_exact=true checks=%o no_RR_or_class_search=true time=%o\n",checks,Cputime(tt);quit;
