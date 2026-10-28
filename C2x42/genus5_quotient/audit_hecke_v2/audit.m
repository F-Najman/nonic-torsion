SetColumns(0);SetQuitOnError(true);SetMemoryLimit(7500*10^6);tt:=Cputime();pr:=5;
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
load "common/native_group_generators.m";
load "common/class_group_direct_RR_checked_v4_fastfb.m";
 k:=GF(pr);AV<v>:=FunctionField(k);bv:=v^2*(v-1);cv:=v*(v-1);EV:=EllipticCurve([AV|1-cv,-bv,-bv,0,0]);QV:=EV![0,0];assert 7*QV eq EV!0 and QV ne EV!0 and 2*QV eq EV![bv,bv*cv];
 hv:=v^2-v+1;sv:=(v^3-3*v+1)/(v*(v-1));tv:=sv-8;dv:=tv^2+13*tv+49;nv:=tv^2+5*tv+1;mv:=tv^4+14*tv^3+63*tv^2+70*tv-7;
 b2:=(1-cv)^2-4*bv;b4:=-(1-cv)*bv;b6:=bv^2;c4:=b2^2-24*b4;c6:=-b2^3+36*b2*b4-216*b6;
 assert dv eq hv^3/(v^2*(v-1)^2) and c4/hv^4 eq nv/dv and c6/hv^6 eq -mv/dv^2 and Discriminant(EV)/hv^12 eq tv/dv^4;
 vp:=(v-1)/v;assert vp^2-vp+1 eq hv/v^2;
 printf "C242_GENUS5_TATE7_DESCENT_IDENTITIES_PASS p=%o cpu=%o\n",pr,Cputime(tt);
 A<u>:=FunctionField(k);t:=u^2;D:=t^2+13*t+49;N:=t^2+5*t+1;M:=t^4+14*t^3+63*t^2+70*t-7;aa:=-N/(48*D);bb:=M/(864*D^2);E:=EllipticCurve([A|0,0,0,aa,bb]);assert Discriminant(E) eq u^2/D^4;
 R<X>:=PolynomialRing(A);psi:=R!DivisionPolynomial(E,3);assert Degree(psi) eq 4;fp:=R![Coefficient(psi,i)/LeadingCoefficient(psi):i in [0..4]];assert IsMonic(fp) and IsIrreducible(fp);K<x3>:=FunctionField(fp);S<Y>:=PolynomialRing(K);fy:=Y^2-(x3^3+(K!aa)*x3+K!bb);assert IsIrreducible(fy);L<y3>:=FunctionField(fy);F:=RationalExtensionRepresentation(L);
 uu:=F!(L!(K!u));xx:=F!(L!x3);yy:=F!y3;assert yy^2 eq xx^3+(F!(L!(K!aa)))*xx+F!(L!(K!bb));assert Genus(F) eq 5;
 jf:=F!(L!(K!jInvariant(E)));cs:=Poles(jf);assert &+[Degree(c):c in cs] eq 8;
 rho:=F!(L!(K!(M/(12*N*D))));hh:=yy/(xx-rho);records:=[* *];branches:=AssociativeArray();for i in [0..3] do branches[i]:=[];end for;
 for j in [1..#cs] do
  P:=cs[j];vu:=Valuation(uu,P);assert vu ne 0 and Valuation(hh,P) eq 0;h2:=Evaluate(hh^2,P);kp:=Parent(h2);
  if vu gt 0 then
   if h2 eq kp!(-1/28) then branch:=1;else assert h2 eq kp!(1/84);branch:=0;end if;
  else
   if h2 eq kp!(1/4) then branch:=2;else assert h2 eq kp!(-1/12);branch:=3;end if;
  end if;
  Append(~branches[branch],j);Append(~records,[*j,Degree(P),vu,Valuation(jf,P),branch,[Integers()!z:z in Coefficients(DefiningPolynomial(kp,k))],[Integers()!z:z in Eltseq(h2,k)]*]);
 end for;assert &and[&+[Degree(cs[j]):j in branches[i]] eq 2:i in [0..3]];assert #branches[2] eq 2 and &and[Degree(cs[j]) eq 1:j in branches[2]];
 base:=Divisor(cs[branches[2][1]]);invs:=ClassGroupAbelianInvariants(F:ReductionDivisor:=base,Proof:=true);MDProvenClassGroupGenerators(F);
 G,_,phi,wit:=MDClassGroupDirectRRCheckedV4FastFB(F,"C2x42/logs/genus5_quotient/audit_hecke_v2/RR.log");hc:=[phi(Divisor(P)-Degree(P)*base):P in cs];hf:=sub<G|hc>;
 hq:=[phi(Divisor(cs[branches[2][2]])-base)] cat [phi(&+[Divisor(cs[j]):j in branches[i]]-2*base):i in [0,1,3]];hqg:=sub<G|hq>;ord:=&*[a:a in Invariants(G)|a ne 0];
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

zcoords:=[uu,xx,yy];c0:=cs[branches[2][1]];assert base eq Divisor(c0);gamma:=[1..#cs];delta:=[1..#cs];
assert #branches[0] eq 2;gamma[branches[0][1]]:=branches[0][2];gamma[branches[0][2]]:=branches[0][1];
for bb in [0,2] do delta[branches[bb][1]]:=branches[bb][2];delta[branches[bb][2]]:=branches[bb][1];end for;
acs0:=[hc[gamma[i]]+11*hc[delta[gamma[i]]]-hc[i]-11*hc[delta[i]]:i in [1..#cs]];baseA:=acs0[branches[2][1]];acs:=[acs0[i]-Degree(cs[i])*baseA:i in [1..#cs]];assert baseA eq G!0;
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
selected:=[];gp:=[];Hg:=hf;
for d in [1..4] do
 for z in Places(F,d) do
  if z in cs or Valuation(jf,z) ne 0 or &or[Valuation(f,z) lt 0:f in zcoords] then continue;end if;
  av:=phi(Divisor(z)-d*base);if av in Hg then continue;end if;
  Append(~selected,z);Append(~gp,av);Hg:=sub<G|hc cat gp>;
  printf "C242_G5_NEW_GENERATOR d=%o order=%o span=%o cpu=%o\n",d,Order(av),#Hg,Cputime(tt);
  if #Hg eq ord then break;end if;
 end for;if #Hg eq ord then break;end if;
end for;assert #Hg eq ord;agp:=[];records:=[* *];proofdata:=[* *];
for i in [1..#selected] do
 z:=selected[i];d:=Degree(z);KD:=GF(5^d);values:=MDM2ToCommonField([Evaluate(f,z):f in zcoords],KD);printf "C242_G5_LIFT_START i=%o d=%o cpu=%o\n",i,d,Cputime(tt);EL,PL,QL:=Lift42(values);print "C242_G5_LIFT_PASS",i,Cputime(tt);assert InverseThree(zcoords,Project42(EL,PL,QL)) eq z;
 allimages:=IsogeniesProof(EL,[PL,QL],11);ipdata:=[* *];assert &+[t[3]:t in allimages] eq 12;Dz:=DivisorGroup(F)!0;imdata:=[* *];
 for q in allimages do
  EG,pts,e,gg,Esource,pointsource,embedK:=Explode(q);vv:=Project42(EG,pts[1],pts[2]);Append(~ipdata,[*e,CurveData(Esource,pointsource[1],pointsource[2]),[ElementData(a):a in Coefficients(gg)],CurveData(EG,pts[1],pts[2]),ElementData(embedK),[ElementData(a):a in vv]*]);zz:=InverseThree(zcoords,vv);numer:=3*d*e;assert numer mod Degree(zz) eq 0;mult:=numer div Degree(zz);Dz+:=mult*Divisor(zz);Append(~imdata,[*Degree(zz),mult,Sprint(zz)*]);
 end for;assert Degree(Dz) eq 36*d;supp,mm:=Support(Dz);assert &and[m mod 3 eq 0:m in mm];Dz:=&+[(mm[j] div 3)*Divisor(supp[j]):j in [1..#supp]];assert Degree(Dz) eq 12*d;dz:=InverseThree(zcoords,[values[1],values[2],-values[3]]);assert Degree(dz) eq d;
 Az:=Dz-Divisor(z)-11*Divisor(dz);assert Degree(Az) eq 0;az:=phi(Az)-d*baseA;Append(~agp,az);Append(~records,[*d,Eltseq(gp[i]),imdata,Eltseq(az)*]);
 Append(~proofdata,[*CurveData(EL,PL,QL),[ElementData(a):a in values],ipdata*]);
 printf "C242_G5_T11_GENERATOR i=%o d=%o class=%o A11=%o cpu=%o\n",i,d,Order(gp[i]),Order(az),Cputime(tt);
end for;
FF:=FreeAbelianGroup(#hc+#gp);hin:=hom<FF->G|hc cat gp>;hout:=hom<FF->G|acs cat agp>;assert #Image(hin) eq ord;assert &and[hout(k) eq G!0:k in Generators(Kernel(hin))];
KR:=sub<G|[hin(k):k in Generators(Kernel(hout))]>;assert #KR*#Image(hout) eq ord;assert &and[a in KR:a in hq];
out:=[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(a):a in hc],[Eltseq(a):a in acs],records,[Eltseq(G!a):a in Generators(KR)],Invariants(KR),Invariants(Image(hout)) *];
Write("C2x42/data/genus5_quotient/audit_hecke_v2/arithmetic.txt",Sprint(out):Overwrite:=true);Write("C2x42/data/genus5_quotient/audit_hecke_v2/explicit_hecke_data.txt",Sprint(proofdata):Overwrite:=true);fprintf wit,"]\n";Flush(wit);
QK,qk:=quo<G|KR>;Write("C2x42/data/genus5_quotient/audit_hecke_v2/quotient_coordinates.txt",Sprint([*Invariants(QK),[Eltseq(r):r in Rows(RelationMatrix(QK))],[Eltseq(qk(G.i)):i in [1..Ngens(G)]]*]):Overwrite:=true);
print "C242_GENUS5_HECKE11_AUDIT_EXPORT_DONE",Invariants(KR),Invariants(Image(hout)),"complete=true",Cputime(tt);quit;
