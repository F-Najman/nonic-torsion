SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "P37/data/input_records.m";
function XYKey(xx,yy)
 d:=LCM(Degree(MinimalPolynomial(xx,GF(2))),Degree(MinimalPolynomial(yy,GF(2))));
 found:=false;
 for bits in [[1,0],[0,1],[1,1]] do z:=bits[1]*xx+bits[2]*yy;mp:=MinimalPolynomial(z,GF(2));if Degree(mp) eq d then bb:=bits;found:=true;break;end if;end for;assert found;
 B:=Matrix(GF(2),[Eltseq(z^i,GF(2)):i in [0..d-1]]);assert Rank(B) eq d;
 cs:=[Eltseq(Solution(B,Vector(GF(2),Eltseq(v,GF(2))))):v in [xx,yy]];
 return [d]cat bb cat[Integers()!c:c in Eltseq(mp)]cat[Integers()!c:c in &cat cs];
end function;
rec:=Inputs[1];assert rec[1] eq 6 and rec[2] eq 1;
PF:=PolynomialRing(GF(2));F<z>:=ext<GF(2)|PF!rec[3]>;
EE:=EllipticCurve([F!PF!v:v in rec[4]]);PN:=EE![F!PF!v:v in rec[5]];
target:={XYKey(xy[1],xy[2]) where xy:=TateNormalForm_xy(EE,k*PN):k in [1..36]};assert #target eq 3;
E:=EllipticCurve("1225h1");fs:=[z[1]:z in Factorization(DivisionPolynomial(E,37))|Degree(z[1]) eq 6];assert #fs eq 3;
E2,phi:=IsogenyFromKernel(E,&*fs);assert Degree(phi) eq 37 and jInvariant(E2) eq -162677523113838677;
seen:={};polys:={};
for f in fs do
 K<x>:=NumberField(f);R<t>:=PolynomialRing(K);a1,a2,a3,a4,a6:=Explode(aInvariants(E));
 L<y>:=ext<K|t^2+(a1*x+a3)*t-(x^3+a2*x^2+a4*x+a6)>;
 EL:=BaseChange(E,L);Q:=EL![L!x,y];assert 37*Q eq EL!0 and Q ne EL!0;
 xy:=[K!v:v in TateNormalForm_xy(EL,Q)];
 f2:=PF![GF(2)!Numerator(c)/GF(2)!Denominator(c):c in Eltseq(f)];assert Degree(f2) eq 6 and IsIrreducible(f2);Include(~polys,f2);
 FF<zz>:=ext<GF(2)|f2>;
 redxy:=[];
 for v in xy do
  vv:=Eltseq(v);assert &and[IsOdd(Denominator(c)):c in vv];
  Append(~redxy,Evaluate(PF![GF(2)!Numerator(c)/GF(2)!Denominator(c):c in vv],zz));
 end for;
 key:=XYKey(redxy[1],redxy[2]);assert key[1] eq 6;Include(~seen,key);
 printf "NONIC37_SEXTIC_LIFT factor_mod2=%o residue_degree=6 global_point_degree=6 time=%o\n",f2,Cputime(tt);
end for;
assert #polys eq 3 and seen eq target;
print "NONIC37_SEXTIC_LIFTS_COMPLETE all_three_residue_classes_lift=true rational37isogeny_verified=true";quit;
