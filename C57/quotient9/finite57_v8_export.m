SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
load "common/native_group_generators.m";
load "common/class_group_direct_RR_checked_v4_fastfb.m";
load "C57/data/quotient9/finite_j_model_data.m";
load "C57/data/quotient9/diamond5_data.m";
K<t>:=FunctionField(GF(2));PF:=PolynomialRing(GF(2));RY<y>:=PolynomialRing(K);
function DecK(c)return (K!PF!c[1])/(K!PF!c[2]);end function;
f:=RY![DecK(c):c in AFPolynomialData];AF<g>:=FunctionField(f);assert Degree(AF) eq AF_degree and Genus(AF) eq 9;
function DecAF(cs)return &+[AF!DecK(cs[i])*g^(i-1):i in [1..#cs]];end function;
Z:=[DecAF(cs):cs in CanonicalCoordinateData];jAF:=DecAF(JFunctionData);cs:=Poles(jAF);assert #cs eq 6 and Degree(cs[1]) eq 1;
function Norm(zz) i:=Index([v ne 0:v in zz],true);assert i ne 0;return [v/zz[i]:v in zz];end function;
function PC(pl)
 vs:=[Valuation(z,pl):z in Z];ind:=Index(vs,Min(vs));return Norm([Evaluate(z/Z[ind],pl):z in Z]);
end function;
cpoints:=[*PC(c):c in cs*];sig:=[];
for i in [1..#cs]do
 z:=cpoints[i];sz:=Norm([&+[Diamond5Matrix[a,b]*z[b]:b in [1..9]]:a in [1..9]]);hits:=[];
 for k in [1..#cs]do
  if Degree(cs[k]) ne Degree(cs[i])then continue;end if;
  kp:=Parent(z[1]);ok,emb:=IsIsomorphic(Parent(cpoints[k][1]),kp);assert ok;zz:=[emb(v):v in cpoints[k]];
  if exists{r:r in [0..Degree(cs[i])-1]|sz eq [v^(2^r):v in zz]}then Append(~hits,k);end if;
 end for;assert #hits eq 1;Append(~sig,hits[1]);
end for;assert [sig[sig[i]]:i in [1..#cs]] eq [1..#cs];
printf "NONIC57_CUSP_DIAMOND_EXACT permutation=%o time=%o\n",sig,Cputime(tt);
base:=Divisor(cs[1]);invs:=ClassGroupAbelianInvariants(AF:ReductionDivisor:=base,Proof:=true);MDProvenClassGroupGenerators(AF);
G,_,phi,wit:=MDClassGroupDirectRRCheckedV4FastFB(AF,"C57/logs/quotient9/finite_RR.log");
hc:=[phi(Divisor(c)-Degree(c)*base):c in cs];alphac:=[hc[sig[i]]-hc[i]-Degree(cs[i])*(hc[sig[1]]-hc[1]):i in [1..#cs]];H:=sub<G|alphac>;Jc:=sub<G|hc>;ord:=&*[v:v in Invariants(G)|v ne 0];assert ord eq 22050;
widths:=[-Valuation(jAF,c):c in cs];sf:=[1..#cs];ww:=[i:i in [1..#cs]|widths[i] eq 57];assert #ww eq 2;sf[ww[1]]:=ww[2];sf[ww[2]]:=ww[1];
Ffree:=FreeAbelianGroup(#cs);beta:=hom<Ffree->H|[H!v:v in alphac]>;
rawA5:=[];
for i in [1..#cs]do v:=Ffree.(sig[sf[i]])+5*Ffree.(sf[i])-Ffree.(sig[i])-5*Ffree.i;Append(~rawA5,v);end for;
a5h:=[rawA5[i]-Degree(cs[i])*rawA5[1]:i in [1..#cs]];
alphaA5:=[&+[G|Eltseq(a5h[i])[j]*alphac[j]:j in [1..#cs]]:i in [1..#cs]];
for rel in Generators(Kernel(beta))do assert &+[G|Eltseq(Ffree!rel)[j]*alphaA5[j]:j in [1..#cs]] eq G!0;end for;
function ImageA5(v)
 r:=Eltseq(v@@beta);return &+[G|r[j]*alphaA5[j]:j in [1..#cs]];
end function;
aop:=hom<H->G|[ImageA5(H.i):i in [1..NumberOfGenerators(H)]]>;HK:=Kernel(aop);assert #HK eq 7;
Write("C57/data/quotient9/finite_group_A5_dual_verified_data.txt",Sprint([*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Degree(c):c in cs],widths,sig,sf,[Eltseq(v):v in hc],[Eltseq(v):v in alphac],Invariants(Jc),Invariants(H),[Eltseq(G!v):v in Generators(HK)],[Eltseq(v):v in alphaA5]*]):Overwrite:=true);
jsets:=AssociativeArray();
load "C57/data/census/N57_f6_records.m";assert N eq 57 and FieldDegree eq 6;
FF:=ext<GF(2)|PolynomialRing(GF(2))!FieldPolynomial>;jsets[6]:={};
for rr in Records do jj:=Evaluate(PolynomialRing(GF(2))!rr[1],FF.1);Include(~jsets[6],Eltseq(MinimalPolynomial(jj,GF(2))));end for;
load "C57/data/census/N57_f7_records.m";assert N eq 57 and FieldDegree eq 7;
FF:=ext<GF(2)|PolynomialRing(GF(2))!FieldPolynomial>;jsets[7]:={};
for rr in Records do jj:=Evaluate(PolynomialRing(GF(2))!rr[1],FF.1);Include(~jsets[7],Eltseq(MinimalPolynomial(jj,GF(2))));end for;
load "C57/data/census/N57_f8_records.m";assert N eq 57 and FieldDegree eq 8;
FF:=ext<GF(2)|PolynomialRing(GF(2))!FieldPolynomial>;jsets[8]:={};
for rr in Records do jj:=Evaluate(PolynomialRing(GF(2))!rr[1],FF.1);Include(~jsets[8],Eltseq(MinimalPolynomial(jj,GF(2))));end for;
load "C57/data/census/N57_f9_records.m";assert N eq 57 and FieldDegree eq 9;
FF:=ext<GF(2)|PolynomialRing(GF(2))!FieldPolynomial>;jsets[9]:={};
for rr in Records do jj:=Evaluate(PolynomialRing(GF(2))!rr[1],FF.1);Include(~jsets[9],Eltseq(MinimalPolynomial(jj,GF(2))));end for;
function Key(z,dd)
 zz:=Norm(z);return Min([[Integers()!b:b in Eltseq(v^(2^a),GF(2)),v in zz]:a in [0..dd-1]]);
end function;
rawrecords:=[* *];kept:=[* *];allplaces:=[* *];
for ee in [1..9]do
 pls:=Places(AF,ee);FF:=GF(2^ee);pv:=[* *];keys:=AssociativeArray();
 for k in [1..#pls]do
  pl:=pls[k];z:=PC(pl);ok,emb:=IsIsomorphic(Parent(z[1]),FF);assert ok;z:=[emb(v):v in z];Append(~pv,z);Append(~allplaces,[*ee,k,Eltseq(DefiningPolynomial(FF,GF(2))),[[Integers()!v:v in Eltseq(w,GF(2))]:w in z],Valuation(jAF,pl) lt 0*]);ky:=Key(z,ee);assert not IsDefined(keys,ky);keys[ky]:=k;
 end for;
 for kk in [1..#pls]do
  pl:=pls[kk];if Valuation(jAF,pl) lt 0 then continue;end if;
  jj:=Evaluate(jAF,pl);jp:=Eltseq(MinimalPolynomial(jj,GF(2)));ds:=[dd:dd in [6..9]|dd mod ee eq 0 and jp in jsets[dd]];
  if #ds eq 0 then continue;end if;
  z:=pv[kk];sz:=Norm([&+[Diamond5Matrix[a,b]*z[b]:b in [1..9]]:a in [1..9]]);sk:=Key(sz,ee);assert IsDefined(keys,sk);sq:=keys[sk];
  cls:=phi(Divisor(pls[sq])-Divisor(pl)-ee*(Divisor(cs[sig[1]])-base));assert cls in H;
  for dd in ds do ww:=dd div ee;yes:=ww*cls in HK;rec:=[*dd,ee,ww,jp,kk,sq,Eltseq(cls),yes,[Integers()!v:v in Eltseq(DefiningPolynomial(FF,GF(2)))],[[Integers()!v:v in Eltseq(w,GF(2))]:w in z]*];Append(~rawrecords,rec);if yes then Append(~kept,rec);end if;end for;
 end for;
 printf "NONIC57_TARGET_DEGREE_TESTED d=%o places=%o tested=%o kept=%o time=%o\n",ee,#pls,#rawrecords,#kept,Cputime(tt);
end for;
Write("C57/data/quotient9/noncuspidal_filter_v2_data.txt",Sprint([*rawrecords,kept,allplaces*]):Overwrite:=true);
printf "NONIC57_QUOTIENT9_NC_FILTER tested=%o kept=%o complete=true time=%o\n",#rawrecords,#kept,Cputime(tt);
SS:=CuspidalSubspace(ModularSymbolsH(57,[4],2,1));TT:=HeckeOperator(SS,2);DD:=(TT^2-HeckeOperator(SS,4))/2;assert DD eq IdentityMatrix(Rationals(),9);
RR<zz>:=PolynomialRing(Rationals());QQ:=FieldOfFractions(RR);cp:=CharacteristicPolynomial(TT);FP:=RR!(zz^9*Evaluate(cp,(QQ!zz^2+2)/(QQ!zz)));assert Degree(FP) eq 18;
Write("C57/data/quotient9/quotient_frobenius2.txt",Sprint([Integers()!v:v in Eltseq(FP)]):Overwrite:=true);
fprintf wit,"]\n";Flush(wit);
printf "NONIC57_FINITE_ALPHA_A5_KERNEL G=%o cuspgroup=%o alpha_cusps=%o A5kernel=%o widths=%o frob5=%o time=%o\n",Invariants(G),Invariants(Jc),Invariants(H),Invariants(HK),widths,sf,Cputime(tt);print "NONIC57_QUOTIENT9_NC_EXPORT_COMPLETE";quit;
