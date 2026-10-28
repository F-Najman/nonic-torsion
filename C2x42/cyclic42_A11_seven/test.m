SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");Attach("common/fast_hecke_cyclic_commonzero.m");Attach("common/class_rational_reduction.m");
load "common/pair_af_place_linear.m";
Y:=MDX1(42,GF(5));assert Genus(Y) eq 25;C:=Curve(Y);cs:=Cusps(Y);os:=CuspOrbitsQ(Y);c0:=Support([o:o in os|Degree(o) eq 1][1])[1];AF,am:=AlgorithmicFunctionField(FunctionField(C));base:=FunctionFieldDivisor(Divisor(c0));coords:=[am(FunctionField(C).i):i in [1,2]];g:=25;expo:=131040;
csdata:=[* *];split:=[];for c in cs do
 d:=Degree(c);if d gt 5 then continue;end if;oo:=[o:o in os|c in Support(o)];assert #oo eq 1;od:=Degree(oo[1]);ww:=-Valuation(jInvariantMap(Y),c);Append(~csdata,[*d,od,ww*]);if d ne od then assert d eq 3 and od eq 6 and ww in [1,2];Append(~split,c);end if;
end for;assert #split eq 4;Write("C2x42/data/cyclic42_A11_seven/cusp_audit.txt",Sprint(csdata):Overwrite:=true);print "C242_SCALAR_CUSP_AUDIT_PASS",csdata,Cputime(tt);
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
OF:=MaximalOrderFinite(AF);OI:=MaximalOrderInfinite(AF);bf:=[AF!b:b in Basis(OF)];bi:=[AF!b:b in Basis(OI)];Write("C2x42/data/cyclic42_A11_seven/model_bases.txt",Sprint([*Degree(AF),[CoeffData(c):c in Coefficients(DefiningPolynomial(AF))],[FData(c):c in coords],[FData(b):b in bf],[FData(b):b in bi]*]):Overwrite:=true);print "C242_A11_SEVEN_MODEL_BASES_READY",Degree(AF),Cputime(tt);
input:=eval Read("C2x42/data/filter_quotient14_linear_arithmetic.txt");fields:=input[2];own:=AssociativeArray();sm:=AssociativeArray();for d in [4..6] do own[d]:=[r:r in input[3]|r[1] eq d];end for;for d in [4,5] do sm[d]:=[r:r in input[4]|r[1] eq d];end for;
R0:=PolynomialRing(GF(5));fs:=AssociativeArray();for fd in fields do fs[fd[1]]:=ext<GF(5)|R0!fd[2]>;end for;
cache:=AssociativeArray();pointProof:=[* *];cases:=[*[1,4,5],[1,6,39],[2,4,3,4,55],[2,4,3,4,65],[2,4,3,5,260],[2,4,3,5,461],[2,4,3,5,589]*];records:=[* *];left:=[* *];
procedure A11Point(~answer,~cache,~pointProof,row)
 d:=row[1];ky:=row[3];kk:=[d] cat ky;if IsDefined(cache,kk) then answer:=cache[kk];return;end if;
 K:=fs[d];z:=K.1;b:=Evaluate(R0!ky[1..d],z);c:=Evaluate(R0!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Order(Q) eq 42 and P ne 21*Q;
 y:=C48PairAFPlaceLinear(coords,TateNormalForm_xy(E,Q));yd:=C48PairAFPlaceLinear(coords,TateNormalForm_xy(E,11*Q));de:=Degree(y);assert de eq Degree(yd) and d mod de eq 0;weight:=d div de;isog:=IsogeniesProof(E,[P,Q],11);assert &+[q[3]:q in isog] eq 12;Dz:=DivisorGroup(AF)!0;idata:=[* *];
 for q in isog do
  EG,pts,e,gg,Esource,pointsource,embedK:=Explode(q);vv:=TateNormalForm_xy(EG,pts[2]);yy:=C48PairAFPlaceLinear(coords,vv);dd:=Degree(yy);assert d*e mod dd eq 0;mult:=d*e div dd;Dz+:=mult*Divisor(yy);Append(~idata,[*e,CurveData(Esource,pointsource[1],pointsource[2]),[ElementData(a):a in Coefficients(gg)],CurveData(EG,pts[1],pts[2]),ElementData(embedK),[ElementData(a):a in vv],dd,mult*]);
 end for;assert Degree(Dz) eq 12*d;answer:=Dz-weight*Divisor(y)-11*weight*Divisor(yd);assert Degree(answer) eq 0;cache[kk]:=answer;Append(~pointProof,[*kk,CurveData(E,P,Q),de,weight,idata*]);Write("C2x42/data/cyclic42_A11_seven/explicit_hecke_data.txt",Sprint(pointProof):Overwrite:=true);print "C242_A11_SEVEN_POINT_READY",kk[1],#pointProof,Cputime(tt);
end procedure;
for j in [1..#cases] do
 desc:=cases[j];rows:=[*own[desc[2]][desc[3]]*];if desc[1] eq 2 then Append(~rows,sm[desc[4]][desc[5]]);end if;Adiv:=DivisorGroup(AF)!0;
 for row in rows do a:=DivisorGroup(AF)!0;A11Point(~a,~cache,~pointProof,row);Adiv+:=a;end for;assert Degree(Adiv) eq 0;print "C242_A11_SEVEN_CASE_RR_START",j,desc,Cputime(tt);bbasis:=Basis(Adiv:Reduction:=false,Simplification:="None");assert #bbasis in [0,1];if #bbasis eq 1 then assert Divisor(bbasis[1])+Adiv eq DivisorGroup(AF)!0;Append(~left,desc);end if;
 I:=ideal<OF|1>;J:=ideal<OI|1>;ss,mm:=Support(Adiv);for k in [1..#ss] do if IsFinite(ss[k]) then I*:=Ideal(ss[k])^(-mm[k]);else J*:=Ideal(ss[k])^(-mm[k]);end if;end for;
 Write(Sprintf("C2x42/data/cyclic42_A11_seven/ideal_bases_case%o.txt",j),Sprint([*desc,[FData(b):b in Basis(I)],[FData(b):b in Basis(J)]*]):Overwrite:=true);
 Append(~records,[*j,desc,#bbasis,Sprint(Adiv),[Sprint(f):f in bbasis]*]);Write("C2x42/data/cyclic42_A11_seven/arithmetic.txt",Sprint([*records,left*]):Overwrite:=true);print "C242_A11_SEVEN_CASE_DONE",j,#bbasis,Cputime(tt);
end for;print "C242_A11_SEVEN_DONE remaining",left,"complete=true",Cputime(tt);quit;
