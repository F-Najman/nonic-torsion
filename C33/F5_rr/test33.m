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
out:=Open("C33/data/F5_rr/dimensions33.jsonseq","w");dims:=[];
for i in[1..#CuspRepresentatives]do
 C:=&+[CuspRepresentatives[i][j]*blocks[j]:j in[1..17]];assert Degree(C)eq 9;
 V,mp:=RiemannRochSpace(C);dim:=Dimension(V);Append(~dims,dim);fprintf out,"{\"index\":%o,\"class\":%o,\"degree\":9,\"h0\":%o}\n",i,ClassLabels[i],dim;Flush(out);
 printf "NONIC33_F5_UNKNOWN_CLASS_RR index=%o h0=%o time=%o\n",i,dim,Cputime(tt);
end for;delete out;assert not assigned F`cld;
printf "NONIC33_F5_RR_TESTS_COMPLETE dims=%o allzero=%o no_classgroup=true time=%o\n",dims,&and[v eq 0:v in dims],Cputime(tt);quit;
