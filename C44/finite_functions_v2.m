SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");Attach("common/fast_hecke.m");
load "C44/data/canonical_data.m";
assert AFBaseCoordinate eq 2;PF:=PolynomialRing(GF(3));X:=MDX1(44,GF(3));C:=Curve(X);cf:=DefiningPolynomial(C);AA<xx,yy>:=AffineSpace(GF(3),2);AR:=CoordinateRing(AA);
if Rank(Parent(cf)) eq 3 then fxy:=Evaluate(cf,[AR.1,AR.2,AR!1]);else fxy:=Evaluate(cf,[AR.1,AR.2]);end if;fx:=Derivative(fxy,1);fy:=Derivative(fxy,2);
function Normalize(z)
 k:=Index([w ne 0:w in z],true);assert k gt 0;return [w/z[k]:w in z];
end function;
function CanonicalValues(xy)
 LF:=Parent(xy[1]);tv:=xy[2];gv:=xy[1];gp:=[LF!1];for i in [1..18] do Append(~gp,gp[#gp]*gv);end for;zz:=[];
 for co in CanonicalCoordinateData do val:=LF!0;for j in [1..#co] do val+:=(Evaluate(PF!co[j][1],tv)/Evaluate(PF!co[j][2],tv))*gp[j];end for;Append(~zz,val);end for;return Normalize(zz);
end function;
function IntCode(v)
 z:=0;pw:=1;for c in v do z+:=(Integers()!c)*pw;pw*:=3;end for;return z;
end function;
function ReduceXY(xy)
 LF:=Parent(xy[1]);DD:=Degree(LF);dd:=Min([r:r in Divisors(DD)|&and[w^(3^r) eq w:w in xy]]);
 candidates:=[xy[1],xy[2],xy[1]+xy[2],xy[1]+2*xy[2]];ix:=0;
 for i in [1..4] do mp:=MinimalPolynomial(candidates[i],GF(3));if Degree(mp) eq dd then ix:=i;theta:=candidates[i];break;end if;end for;assert ix gt 0;
 BB:=Matrix(GF(3),[Eltseq(theta^i,GF(3)):i in [0..dd-1]]);assert Nrows(BB) eq dd and Ncols(BB) eq DD and Rank(BB) eq dd;
 cx:=Eltseq(Solution(BB,Vector(GF(3),Eltseq(xy[1],GF(3)))));cy:=Eltseq(Solution(BB,Vector(GF(3),Eltseq(xy[2],GF(3)))));
 assert Evaluate(PF!cx,theta) eq xy[1] and Evaluate(PF!cy,theta) eq xy[2];
 RF<ar>:=ext<GF(3)|mp>;xy:=[Evaluate(PF!cx,ar),Evaluate(PF!cy,ar)];assert Degree(RF) eq dd and Evaluate(fxy,xy) eq 0;
 z:=CanonicalValues(xy);key:=<dd,ix,IntCode(Eltseq(mp)),IntCode(cx),IntCode(cy)>;return xy,z,dd,key;
end function;
function EncodingBasis(LF)
 dd:=Degree(LF);theta:=LF.1;if Degree(MinimalPolynomial(theta,GF(3))) ne dd then theta:=PrimitiveElement(LF);end if;
 BB:=Matrix(GF(3),[Eltseq(theta^i,GF(3)):i in [0..dd-1]]);assert Nrows(BB) eq dd and Ncols(BB) eq dd and Rank(BB) eq dd;return [Integers()!c:c in Eltseq(MinimalPolynomial(theta,GF(3)))],BB^-1;
end function;
function Enc(w,IB)return [Integers()!a:a in Eltseq(Vector(GF(3),Eltseq(w,GF(3)))*IB)];end function;
function JetData(xy,zz,mm)
 LF:=Parent(xy[1]);prec:=Max(mm,2);PS<T>:=PowerSeriesRing(LF,prec);xs:=PS!xy[1];ys:=PS!xy[2];vx:=Evaluate(fx,xy);vy:=Evaluate(fy,xy);
 if vy ne 0 then dir:=1;xs+:=T;for k in [1..prec-1] do er:=Coefficient(Evaluate(fxy,[xs,ys]),k);ys-:=(er/vy)*T^k;end for;
 else assert vx ne 0;dir:=2;ys+:=T;for k in [1..prec-1] do er:=Coefficient(Evaluate(fxy,[xs,ys]),k);xs-:=(er/vx)*T^k;end for;end if;
 assert &and[Coefficient(Evaluate(fxy,[xs,ys]),k) eq 0:k in [0..prec-1]];
 gp:=[PS!1];for i in [1..18] do Append(~gp,gp[#gp]*xs);end for;zs:=[];
 for co in CanonicalCoordinateData do val:=PS!0;for j in [1..#co] do val+:=(Evaluate(PF!co[j][1],ys)/Evaluate(PF!co[j][2],ys))*gp[j];end for;Append(~zs,val);end for;
 ki:=Index([Coefficient(w,0) ne 0:w in zs],true);assert ki gt 0;zs:=[w/zs[ki]:w in zs];assert [Coefficient(w,0):w in zs] eq zz and exists{w:w in zs|Coefficient(w,1) ne 0};
 fp,IB:=EncodingBasis(LF);xycode:=[Enc(w,IB):w in xy];xyjets:=[[Enc(Coefficient(w,k),IB):k in [0..mm-1]]:w in [xs,ys]];zjets:=[[Enc(Coefficient(w,k),IB):k in [0..mm-1]]:w in zs];return fp,xycode,xyjets,zjets,dir;
end function;
