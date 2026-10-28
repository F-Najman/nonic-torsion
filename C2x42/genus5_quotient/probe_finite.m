/* New low-genus quotient probe; no global torsion conclusion is assumed. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(7500*10^6);tt:=Cputime();
load "common/native_group_generators.m";
load "common/class_group_direct_RR_checked_v4_fastfb.m";
function Probe(pr)
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
 G,_,phi,wit:=MDClassGroupDirectRRCheckedV4FastFB(F,Sprintf("C2x42/logs/genus5_quotient/p%o_RR.log",pr));hc:=[phi(Divisor(P)-Degree(P)*base):P in cs];hf:=sub<G|hc>;
 hq:=[phi(Divisor(cs[branches[2][2]])-base)] cat [phi(&+[Divisor(cs[j]):j in branches[i]]-2*base):i in [0,1,3]];hqg:=sub<G|hq>;ord:=&*[a:a in Invariants(G)|a ne 0];
 data:=[*pr,Genus(F),Degree(K),Degree(L),records,Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(a):a in hc],Invariants(hf),[Eltseq(a):a in hq],Invariants(hqg),ord*];
 Write(Sprintf("C2x42/data/genus5_quotient/p%o_arithmetic.txt",pr),Sprint(data):Overwrite:=true);fprintf wit,"]\n";Flush(wit);
 printf "C242_GENUS5_PROVEN_FINITE_GROUP p=%o G=%o Hfinite=%o Hrational_cusps_candidate=%o order=%o cpu=%o\n",pr,Invariants(G),Invariants(hf),Invariants(hqg),ord,Cputime(tt);return data;
end function;
data:=[*Probe(5),Probe(11)*];Write("C2x42/data/genus5_quotient/combined_arithmetic.txt",Sprint(data):Overwrite:=true);print "C242_GENUS5_QUOTIENT_FINITE_PROBE_DONE complete=true",Cputime(tt);quit;
