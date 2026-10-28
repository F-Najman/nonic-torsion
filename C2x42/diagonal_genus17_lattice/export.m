SetColumns(0);SetQuitOnError(true);SetMemoryLimit(7500*10^6);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");Attach("common/fast_diamond_hash_m2.m");Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
function InverseThree(coords,values)
 F:=Parent(coords[1]);K:=Parent(values[1]);k:=PrimeField(K);p:=#k;e:=LCM([Degree(MinimalPolynomial(v,k)):v in values]);found:=false;
 for bits in [[1,0,0],[0,1,0],[0,0,1]] cat [[1,a,b]:a,b in [0..p-1]] cat [[0,1,a]:a in [0..p-1]] do
  gen:=&+[K!bits[i]*values[i]:i in [1..3]];pol:=MinimalPolynomial(gen,k);if Degree(pol) eq e then chosen:=bits;found:=true;break;end if;
 end for;assert found;M:=Matrix(k,[Eltseq(gen^i,k):i in [0..e-1]]);assert Rank(M) eq e;
 cc:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in values];assert &and[&+[K!cc[j][i+1]*gen^i:i in [0..e-1]] eq values[j]:j in [1..3]];
 gg:=&+[F!chosen[i]*coords[i]:i in [1..3]];fs:=[Evaluate(pol,gg)] cat [coords[j]-&+[F!cc[j][i+1]*gg^i:i in [0..e-1]]:j in [1..3]];
 O:=MaximalOrderFinite(F);I:=ideal<O|fs>;fa:=Factorization(I);places:=[Place(t[1]):t in fa|t[2] gt 0];
 if #places eq 0 then Oi:=MaximalOrderInfinite(F);Ii:=ideal<Oi|fs>;fi:=Factorization(Ii);places:=[Place(t[1]):t in fi|t[2] gt 0];end if;
 assert #places eq 1 and Degree(places[1]) eq e;assert &and[f eq 0 or Valuation(f,places[1]) gt 0:f in fs];return places[1];
end function;

function Project42(E,P,Q)
 aa:=aInvariants(E);assert aa[1] eq 0 and aa[3] eq 0;
 S:=6*Q;T:=14*Q;R:=21*Q;Z:=P+R;assert 7*S eq E!0 and 3*T eq E!0;
 xs:=S[1];ys:=S[2];m:=(3*xs^2+2*aa[2]*xs+aa[4])/(2*ys);a2p:=3*xs+aa[2]-m^2;a3p:=2*ys;
 kv:=a3p/a2p;bv:=-a2p^3/a3p^2;cv:=1-2*m*a2p/a3p;v:=bv/cv;assert bv eq v^2*(v-1) and cv eq v*(v-1);H:=v^2-v+1;assert H ne 0;
 u:=4*(P[1]-R[1])*(P[1]-Z[1])*(R[1]-Z[1])/(kv^6*v^4*(v-1)^4);
 xx:=(T[1]+aa[2]/3)/(kv^2*H^2);yy:=T[2]/(kv^3*H^3);
 t:=u^2;D:=t^2+13*t+49;N:=t^2+5*t+1;M:=t^4+14*t^3+63*t^2+70*t-7;
 assert yy^2 eq xx^3-N/(48*D)*xx+M/(864*D^2);assert 7*Q ne E!0 and 42*Q eq E!0;
 return [u,xx,yy];
end function;
function Lift42(values)
 d:=Degree(Parent(values[1]));L:=GF(5^(3*d));u,x,y:=Explode([L!z:z in values]);t:=u^2;D:=t^2+13*t+49;N:=t^2+5*t+1;M:=t^4+14*t^3+63*t^2+70*t-7;
 EL:=EllipticCurve([L|0,0,0,-N/(48*D),M/(864*D^2)]);T:=EL![x,y];assert 3*T eq EL!0 and T ne EL!0;
 R<W>:=PolynomialRing(L);s:=t+8;vr:=Roots(W^3-s*W^2+(s-3)*W+1);assert #vr eq 3;v:=vr[1][1];bv:=v^2*(v-1);cv:=v*(v-1);H:=v^2-v+1;assert H ne 0;
 a2v:=-bv+(1-cv)^2/4;S:=EL![a2v/(3*H^2),-bv/(2*H^3)];assert Order(S) eq 7;
 tw:=Roots(W^3-N/(48*D)*W+M/(864*D^2));assert #tw eq 3;P:=EL![tw[1][1],0];R2:=EL![tw[2][1],0];Q:=R2+2*T+6*S;vv:=Project42(EL,P,Q);
 if vv[1] ne u then P,R2:=Explode([R2,P]);Q:=R2+2*T+6*S;vv:=Project42(EL,P,Q);end if;
 assert vv eq [u,x,y] and Order(Q) eq 42 and P ne 21*Q;
 return EL,P,Q;
end function;

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

k:=GF(5);BU<u>:=FunctionField(k);t:=u^2;D:=t^2+13*t+49;N:=t^2+5*t+1;M:=t^4+14*t^3+63*t^2+70*t-7;RU<W>:=PolynomialRing(BU);fw:=W^3-N/48*W+M*(2*t+13)/(1728*D)-u/(8*D);assert IsIrreducible(fw);K0<w>:=FunctionField(fw);u0:=K0!u;assert Genus(K0) eq 0;ps:=Poles(u0);assert #ps eq 3 and &and[Degree(p) eq 1:p in ps];bbasis:=Basis(Divisor(ps[1]));assert #bbasis eq 2;z0:=[v:v in bbasis|Valuation(v,ps[1]) lt 0][1];assert &+[Degree(p)*(-Valuation(z0,p)):p in Poles(z0)] eq 1;
vs:=[z0^i:i in [0..3]] cat [u0*z0^i:i in [0..2]];vv:=[Eltseq(v):v in vs];den:=LCM(&cat[[Denominator(c):c in row]:row in vv]);deg:=Maximum(&cat[[Degree(Numerator(c*den)):c in row]:row in vv]);mat:=Matrix(k,[[Coefficient(Numerator(row[j]*den),i):j in [1..3],i in [0..deg]]:row in vv]);ker:=Nullspace(mat);assert Dimension(ker) eq 1;rel:=Eltseq(Basis(ker)[1]);RZ<Z>:=PolynomialRing(k);pn:=RZ!rel[1..4];pd:=RZ!rel[5..7];assert Evaluate(pn,z0)+u0*Evaluate(pd,z0) eq 0;
BZ<z>:=FunctionField(k);uz:=-Evaluate(pn,z)/Evaluate(pd,z);assert Maximum(Degree(Numerator(uz)),Degree(Denominator(uz))) eq 3;print "C242_DIAGONAL_GENUS0_PARAMETRIZATION",Numerator(uz),Denominator(uz),Cputime(tt);
function CoeffData(f)
 return [[Integers()!a:a in Coefficients(Numerator(f))],[Integers()!a:a in Coefficients(Denominator(f))]];
end function;
zdata:=[CoeffData(c):c in Eltseq(z0)];udata:=CoeffData(uz);Write("C2x42/data/diagonal_genus17_lattice/parameter_data.txt",Sprint([*zdata,udata*]):Overwrite:=true);
t:=uz^2;D:=t^2+13*t+49;N:=t^2+5*t+1;M:=t^4+14*t^3+63*t^2+70*t-7;aa:=-N/(48*D);bb:=M/(864*D^2);E0:=EllipticCurve([BZ|0,0,0,aa,bb]);RX<X>:=PolynomialRing(BZ);psi:=RX!DivisionPolynomial(E0,3);fp:=RX![Coefficient(psi,i)/LeadingCoefficient(psi):i in [0..4]];assert IsMonic(fp) and IsIrreducible(fp);K1<x3>:=FunctionField(fp);RY<Y>:=PolynomialRing(K1);fy:=Y^2-(x3^3+(K1!aa)*x3+K1!bb);assert IsIrreducible(fy);L<y3>:=FunctionField(fy);F:=RationalExtensionRepresentation(L);zzF:=F!(L!(K1!z));xxF:=F!(L!x3);yyF:=F!y3;coords17:=[zzF,xxF,yyF];assert Genus(F) eq 17;jf:=F!(L!(K1!jInvariant(E0)));cs:=Poles(jf);assert &+[Degree(c):c in cs] eq 24;print "C242_DIAGONAL_GENUS17_MODEL_READY",Genus(F),[Degree(c):c in cs],Cputime(tt);
Write("C2x42/data/diagonal_genus17_lattice/normalized_model_data.txt",Sprint([*[CoeffData(c):c in Coefficients(DefiningPolynomial(F))],[[CoeffData(c):c in Eltseq(f)]:f in coords17]*]):Overwrite:=true);
function EvalParameter(zdata,u,w)
 K:=Parent(u);R:=PolynomialRing(K);return &+[Evaluate(R!zdata[i][1],u)/Evaluate(R!zdata[i][2],u)*w^(i-1):i in [1..#zdata]];
end function;
function Project17(E,P,Q,zdata)
 vals:=Project42(E,P,Q);S:=6*Q;ai:=aInvariants(E);x:=S[1];y:=S[2];m:=(3*x^2+2*ai[2]*x+ai[4])/(2*y);a2p:=3*x+ai[2]-m^2;a3p:=2*y;kv:=a3p/a2p;b:=-a2p^3/a3p^2;c:=1-2*m*a2p/a3p;v0:=b/c;H:=v0^2-v0+1;v1:=(v0-1)/v0;v2:=1/(1-v0);ts:=[P,21*Q,P+21*Q];rr:=[(T[1]+ai[2]/3)/(kv^2*H^2):T in ts];w:=rr[1]*v0+rr[2]*v1+rr[3]*v2;u:=vals[1];t:=u^2;D:=t^2+13*t+49;N:=t^2+5*t+1;M:=t^4+14*t^3+63*t^2+70*t-7;assert w^3-N/48*w+M*(2*t+13)/(1728*D)-u/(8*D) eq 0;
 return [EvalParameter(zdata,u,w),vals[2],vals[3]];
end function;
input:=eval Read("C2x42/data/filter_quotient14_linear_arithmetic.txt");row:=[r:r in input[3]|r[1] eq 9][8];fd:=[a:a in input[2]|a[1] eq 9][1];K:=ext<GF(5)|PolynomialRing(GF(5))!fd[2]>;ky:=row[3];a:=K.1;R0:=PolynomialRing(GF(5));b:=Evaluate(R0!ky[1..9],a);c:=Evaluate(R0!ky[10..18],a);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Order(Q) eq 42 and P ne 21*Q;vals:=Project17(E,P,Q,zdata);assert Project17(E,21*Q,P+37*Q,zdata) eq vals;place:=InverseThree(coords17,vals);assert Degree(place) eq 3;deltaPlace:=InverseThree(coords17,Project17(E,P,11*Q,zdata));assert Degree(deltaPlace) eq 3;
Write("C2x42/data/diagonal_genus17_lattice/source_data.txt",Sprint([*CurveData(E,P,Q),[ElementData(v):v in vals],[ElementData(v):v in Project17(E,P,11*Q,zdata)]*]):Overwrite:=true);
images:=IsogeniesProof(E,[P,Q],11);assert &+[t[3]:t in images] eq 12;Dz:=DivisorGroup(F)!0;proof:=[* *];
for q in images do
 EG,pts,e,gg,Esource,pointsource,embedK:=Explode(q);vv:=Project17(EG,pts[1],pts[2],zdata);zz:=InverseThree(coords17,vv);numer:=9*e;assert numer mod Degree(zz) eq 0;Dz+:=(numer div Degree(zz))*Divisor(zz);Append(~proof,[*e,CurveData(Esource,pointsource[1],pointsource[2]),[ElementData(a):a in Coefficients(gg)],CurveData(EG,pts[1],pts[2]),ElementData(embedK),[ElementData(a):a in vv]*]);
end for;assert Degree(Dz) eq 108;supp,mm:=Support(Dz);assert &and[m mod 3 eq 0:m in mm];Dz:=&+[(mm[j] div 3)*Divisor(supp[j]):j in [1..#supp]];assert Degree(Dz) eq 36;Az:=Dz-Divisor(place)-11*Divisor(deltaPlace);assert Degree(Az) eq 0;print "C242_DIAGONAL_A11_DIVISOR_READY",#supp,Cputime(tt);
Write("C2x42/data/diagonal_genus17_lattice/explicit_hecke_data.txt",Sprint(proof):Overwrite:=true);Write("C2x42/data/diagonal_genus17_lattice/A11_divisor.txt",Sprint(Az):Overwrite:=true);Brr:=Basis(Az:Reduction:=false,Simplification:="None");assert #Brr in [0,1];

lam:=(3*xxF^2+(F!(L!(K1!aa))))/(2*yyF);matlam:=Matrix(BZ,[Eltseq(lam^i):i in [0..7]]);assert Rank(matlam) eq 8;
function LamData(f,matlam)
 cc:=Eltseq(Solution(matlam,Vector(BZ,Eltseq(f))));assert &+[F!cc[i]*lam^(i-1):i in [1..8]] eq f;return [CoeffData(c):c in cc];
end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);Ib:=ideal<OF|1>;sA,mA:=Support(Az);for i in [1..#sA] do assert IsFinite(sA[i]);Ib*:=Ideal(sA[i])^(-mA[i]);end for;
bf:=[F!b:b in Basis(OF)];bi:=[F!b:b in Basis(OI)];bI:=[F!b:b in Basis(Ib)];assert #bf eq 8 and #bi eq 8 and #bI eq 8;
Write("C2x42/data/diagonal_genus17_lattice/lattice_bases.txt",Sprint([*[LamData(f,matlam):f in bf],[LamData(f,matlam):f in bi],[LamData(f,matlam):f in bI]*]):Overwrite:=true);
print "C242_GENUS17_LATTICE_EXPORT_DONE",#Brr,"complete=true",Cputime(tt);quit;
