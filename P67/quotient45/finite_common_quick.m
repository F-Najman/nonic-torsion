// Same exact field already certified genus45 by finish_plane_model_certified actualexit0.
// Avoid repeating the expensive genus computation; native decomposition remains exact.
// Model SHA256 12f4666582d58474a7db9301022e20f9cbc8dacb4f98cf087ed27d7f493f03dc
load "P67/data/quotient45/plane_model_certified_data.m";
load "P67/data/quotient15/isogeny3_target_data.m";
K<t>:=FunctionField(GF(2));PF:=PolynomialRing(GF(2));RY<y>:=PolynomialRing(K);
function DecK(c)return (K!PF!c[1])/(K!PF!c[2]);end function;
f:=RY![DecK(c):c in AFPolynomialData];AF<g>:=FunctionField(f);assert Degree(AF) eq AF_degree;
function DecAF(cs)return &+[AF!DecK(cs[i])*g^(i-1):i in [1..#cs]];end function;
if AFBaseCoordinate eq 1 then uAF:=AF!t;vAF:=g;else uAF:=g;vAF:=AF!t;end if;
Z:=[DecAF(cs):cs in CanonicalCoordinateData];jAF:=DecAF(JFunctionData);
assert #Z eq 45;
function CodeK(w)return <[Integers()!a:a in Eltseq(Numerator(w))],[Integers()!a:a in Eltseq(Denominator(w))]>;end function;
function CodeAF(w)return [CodeK(c):c in Eltseq(w)];end function;
function Normalize(z)
 k:=Index([w ne 0:w in z],true);assert k ne 0;return [w/z[k]:w in z];
end function;
function JKey(z,jv)
 z:=Normalize(z);d:=Degree(MinimalPolynomial(jv,GF(2)));F:=Parent(jv);assert #F eq 2^d;
 bb:=Matrix(GF(2),d,d,[c:c in Eltseq(jv^i,GF(2)),i in [0..d-1]]);assert Rank(bb) eq d;
 cc:=[Eltseq(Solution(bb,Vector(GF(2),Eltseq(w,GF(2))))):w in z];
 return <[Integers()!c:c in Eltseq(MinimalPolynomial(jv,GF(2)))],[[Integers()!a:a in row]:row in cc]>;
end function;
function PlaceCoordinates(pl,ZZ)
 vals:=[Valuation(w,pl):w in ZZ];k:=Index(vals,Min(vals));
 return [Evaluate(w/ZZ[k],pl):w in ZZ];
end function;
function SigmaCoordinates(z)
 return Normalize([&+[Diamond2Matrix[i,k]*z[k]:k in [1..45]]:i in [1..45]]);
end function;
