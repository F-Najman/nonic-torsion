SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");Attach("common/fast_hecke_cyclic_commonzero.m");Attach("common/class_rational_reduction.m");
load "common/pair_af_place_linear.m";
Y:=MDX1(42,GF(5));assert Genus(Y) eq 25;C:=Curve(Y);cs:=Cusps(Y);os:=CuspOrbitsQ(Y);c0:=Support([o:o in os|Degree(o) eq 1][1])[1];AF,am:=AlgorithmicFunctionField(FunctionField(C));base:=FunctionFieldDivisor(Divisor(c0));coords:=[am(FunctionField(C).i):i in [1,2]];g:=25;expo:=131040;
csdata:=[* *];split:=[];for c in cs do
 d:=Degree(c);if d gt 5 then continue;end if;oo:=[o:o in os|c in Support(o)];assert #oo eq 1;od:=Degree(oo[1]);ww:=-Valuation(jInvariantMap(Y),c);Append(~csdata,[*d,od,ww*]);if d ne od then assert d eq 3 and od eq 6 and ww in [1,2];Append(~split,c);end if;
end for;assert #split eq 4;Write("C2x42/data/CM_cyclic42_lattice/cusp_audit.txt",Sprint(csdata):Overwrite:=true);print "C242_SCALAR_CUSP_AUDIT_PASS",csdata,Cputime(tt);
function IsogeniesProof(E, points, q)
    assert IsPrime(q) and IsOdd(q) and q ne Characteristic(BaseRing(E));
    K := BaseRing(E);
    fq := DivisionPolynomial(E, q);
    fac := [t[1] : t in Factorization(fq)];
    assert &+[Degree(h) : h in fac] eq (q^2-1) div 2;   // separable: distinct x-coordinates
    a1, a2, a3, a4, a6 := Explode(aInvariants(E));
    unassigned := [true : h in fac];
    result := [* *];
    total := 0;
    for i in [1..#fac] do
        if not unassigned[i] then continue; end if;
        h := fac[i]; d := Degree(h);
        if d eq 1 then F1 := K; x1 := -Coefficient(h, 0)/Coefficient(h, 1);
        else F1 := ext<K | h>; x1 := F1.1; end if;
        // y-coordinate: solve y^2 + a1 x y + a3 y = x^3 + a2 x^2 + a4 x + a6 over F1 or a quadratic extension
        PY<Y> := PolynomialRing(F1);
        ypol := Y^2 + (F1!a1*x1 + F1!a3)*Y - (x1^3 + F1!a2*x1^2 + F1!a4*x1 + F1!a6);
        rts := Roots(ypol);
        if #rts gt 0 then F2 := F1; y1 := rts[1][1];
        else F2 := ext<F1 | 2>; y1 := Roots(PolynomialRing(F2)!ypol)[1][1]; end if;
        E2 := EllipticCurve([F2!a1, F2!a2, F2!a3, F2!a4, F2!a6]);
        P1 := E2![F2!x1, y1];
        assert q*P1 eq E2!0 and P1 ne E2!0;
        xs := [ F1!((k*P1)[1]) : k in [1..(q-1) div 2] ];          // x-coordinates lie in F1
        PX<X> := PolynomialRing(F1);
        g := &*[X - xk : xk in xs];                                  // kernel polynomial of G = <P1>
        cs := Coefficients(g);
        L := sub< F1 | cs cat [F1!K.1] >;                            // field of definition of G, contains K
        e := Degree(L) div Degree(K);
        // mark the factors belonging to this Frobenius-orbit: divisors of Norm_{L/K}(g)
        qK := #K;
        Ng := PX!1; gg := g;
        for j in [1..e] do
            Ng *:= gg;
            gg := PX![c^qK : c in Coefficients(gg)];
        end for;
        assert gg eq g;   // Frobenius^e fixes g
        NgK := PolynomialRing(K)![K!c : c in Coefficients(Ng)];
        deg_assigned := 0;
        for j in [1..#fac] do
            if unassigned[j] and IsDivisibleBy(NgK, fac[j]) then unassigned[j] := false; deg_assigned +:= Degree(fac[j]); end if;
        end for;
        assert deg_assigned eq e*(q-1) div 2;
        // the isogeny over L
        PL<XL> := PolynomialRing(L);
        gL := PL![L!c : c in cs];
        EL := EllipticCurve([L!(F1!a) : a in aInvariants(E)]);
        EG, phi := IsogenyFromKernel(EL, gL);
        images := [phi(EL![L!(F1!(K!P[1])), L!(F1!(K!P[2]))]) : P in points];
        Append(~result, <EG, images, e, gL, EL, [EL![L!(F1!(K!P[1])), L!(F1!(K!P[2]))]:P in points], L!(F1!K.1)>);
        total +:= e;
    end for;
    assert total eq q+1;
    return result;
end function;function FieldData(K)
 return [Integers()!a:a in Coefficients(DefiningPolynomial(K,GF(5)))];
end function;
function ElementData(v)
 return [Integers()!a:a in Eltseq(v,GF(5))];
end function;
function CurveData(E,P,Q)
 return [*FieldData(BaseRing(E)),[ElementData(a):a in aInvariants(E)],[ElementData(a):a in Eltseq(P)],[ElementData(a):a in Eltseq(Q)]*];
end function;



function CoeffData(f)
 return [[Integers()!a:a in Coefficients(Numerator(f))],[Integers()!a:a in Coefficients(Denominator(f))]];
end function;
function FData(f)
 return [CoeffData(c):c in Eltseq(AF!f)];
end function;
OF:=MaximalOrderFinite(AF);OI:=MaximalOrderInfinite(AF);bf:=[AF!b:b in Basis(OF)];bi:=[AF!b:b in Basis(OI)];Write("C2x42/data/CM_cyclic42_lattice/model_bases.txt",Sprint([*Degree(AF),[CoeffData(c):c in Coefficients(DefiningPolynomial(AF))],[FData(c):c in coords],[FData(b):b in bf],[FData(b):b in bi]*]):Overwrite:=true);print "C242_A11_SEVEN_MODEL_BASES_READY",Degree(AF),Cputime(tt);

function GetCurve(data)
 pol,ai,pp,qq:=Explode(data);K:=ext<GF(5)|PolynomialRing(GF(5))!pol>;R:=PolynomialRing(GF(5));z:=K.1;E:=EllipticCurve([Evaluate(R!a,z):a in ai]);P:=E![Evaluate(R!a,z):a in pp];Q:=E![Evaluate(R!a,z):a in qq];return K,E,P,Q;
end function;
proof:=eval Read("C2x42/data/CM_cyclic42_A11/explicit_hecke_data.txt");out:=[* *];
for rec in proof do
 ii,src,images:=Explode(rec);K,E,P,Q:=GetCurve(src);assert Degree(K) eq 6 and Order(Q) eq 42;y:=C48PairAFPlaceLinear(coords,TateNormalForm_xy(E,Q));yd:=C48PairAFPlaceLinear(coords,TateNormalForm_xy(E,11*Q));d:=Degree(y);assert d eq Degree(yd) and 6 mod d eq 0;w:=6 div d;Dz:=DivisorGroup(AF)!0;
 for im in images do
  e,srcL,gdat,target,embed,vv0,dd,mult:=Explode(im);L,EG,PG,QG:=GetCurve(target);R:=PolynomialRing(GF(5));zz:=L.1;vv:=[Evaluate(R!a,zz):a in vv0];assert vv eq TateNormalForm_xy(EG,QG);yy:=C48PairAFPlaceLinear(coords,vv);assert Degree(yy) eq dd and mult*dd eq 6*e;Dz+:=mult*Divisor(yy);
 end for;assert Degree(Dz) eq 72;Adiv:=Dz-w*Divisor(y)-11*w*Divisor(yd);assert Degree(Adiv) eq 0 and #Basis(Adiv:Reduction:=false,Simplification:="None") eq 0;I:=ideal<OF|1>;J:=ideal<OI|1>;ss,mm:=Support(Adiv);for k in [1..#ss] do if IsFinite(ss[k]) then I*:=Ideal(ss[k])^(-mm[k]);else J*:=Ideal(ss[k])^(-mm[k]);end if;end for;
 Write(Sprintf("C2x42/data/CM_cyclic42_lattice/ideal_bases_case%o.txt",ii),Sprint([*ii,[FData(b):b in Basis(I)],[FData(b):b in Basis(J)]*]):Overwrite:=true);Append(~out,[*ii,d,w,0*]);print "C242_CM_LATTICE_EXPORTED",ii,Cputime(tt);
end for;Write("C2x42/data/CM_cyclic42_lattice/arithmetic.txt",Sprint(out):Overwrite:=true);print "C242_CM_LATTICE_EXPORT_DONE complete=true",Cputime(tt);quit;
