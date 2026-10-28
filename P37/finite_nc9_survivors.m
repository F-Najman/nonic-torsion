SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");Attach("common/fast_hecke.m");
tt:=Cputime();N:=37;
function XYKey(E,P)
 xx,yy:=Explode(TateNormalForm_xy(E,P));
 d:=LCM(Degree(MinimalPolynomial(xx,GF(2))),Degree(MinimalPolynomial(yy,GF(2))));
 found:=false;
 for bits in [[1,0],[0,1],[1,1]] do
  z:=bits[1]*xx+bits[2]*yy;mp:=MinimalPolynomial(z,GF(2));
  if Degree(mp) eq d then bb:=bits;found:=true;break;end if;
 end for;assert found;
 B:=Matrix(GF(2),[Eltseq(z^i,GF(2)):i in [0..d-1]]);assert Rank(B) eq d;
 cs:=[Eltseq(Solution(B,Vector(GF(2),Eltseq(v,GF(2))))):v in [xx,yy]];
 return [d] cat bb cat [Integers()!c:c in Eltseq(mp)] cat [Integers()!c:c in &cat cs];
end function;


load "P37/data/input_records.m";
count:=0;total:=0;
for rec in Inputs do
 if rec[1] ne 9 or rec[2] notin [3,4] then continue;end if;
 PF:=PolynomialRing(GF(2));F<z>:=ext<GF(2)|PF!rec[3]>;E:=EllipticCurve([F!PF!v:v in rec[4]]);P:=E![F!PF!v:v in rec[5]];
 assert XYKey(E,4*P) eq XYKey(E,P);orb:={XYKey(E,k*P):k in [1..36]};assert #orb in [1,2];count+:=1;total+:=#orb;
 printf "NONIC37_NC9_SURVIVOR_DELTA4 index=%o orbit_size=%o fixed=true time=%o\n",rec[2],#orb,Cputime(tt);
end for;assert count eq 2 and total eq 3;
print "NONIC37_NC9_SURVIVORS_ALL_DELTA4_FIXED";quit;
