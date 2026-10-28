SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
load "common/native_group_generators.m";
load "common/class_group_direct_RR_checked_v4_fastfb.m";
load "C65/data/quotient9/finite_j_model_data.m";
load "C65/data/quotient9/diamond3_data.m";
K<t>:=FunctionField(GF(2));PF:=PolynomialRing(GF(2));RY<y>:=PolynomialRing(K);
function DecK(c)return (K!PF!c[1])/(K!PF!c[2]);end function;
f:=RY![DecK(c):c in AFPolynomialData];AF<g>:=FunctionField(f);assert Degree(AF) eq AF_degree and Genus(AF) eq 9;
function DecAF(cs)return &+[AF!DecK(cs[i])*g^(i-1):i in [1..#cs]];end function;
Z:=[DecAF(cs):cs in CanonicalCoordinateData];jAF:=DecAF(JFunctionData);cs:=Poles(jAF);assert #cs eq 6;basenum:=[i:i in[1..#cs]|Degree(cs[i]) eq 1 and Valuation(jAF,cs[i]) eq -1][1];
function Norm(zz) i:=Index([v ne 0:v in zz],true);assert i ne 0;return [v/zz[i]:v in zz];end function;
function PC(pl)
 vs:=[Valuation(z,pl):z in Z];ind:=Index(vs,Min(vs));return Norm([Evaluate(z/Z[ind],pl):z in Z]);
end function;
cpoints:=[*PC(c):c in cs*];sig:=[];
for i in [1..#cs]do
 z:=cpoints[i];sz:=Norm([&+[Diamond3Matrix[a,b]*z[b]:b in [1..9]]:a in [1..9]]);hits:=[];
 for k in [1..#cs]do
  if Degree(cs[k]) ne Degree(cs[i])then continue;end if;
  kp:=Parent(z[1]);ok,emb:=IsIsomorphic(Parent(cpoints[k][1]),kp);assert ok;zz:=[emb(v):v in cpoints[k]];
  if exists{r:r in [0..Degree(cs[i])-1]|sz eq [v^(2^r):v in zz]}then Append(~hits,k);end if;
 end for;assert #hits eq 1;Append(~sig,hits[1]);
end for;assert [sig[sig[i]]:i in [1..#cs]] eq [1..#cs];
printf "NONIC65_CUSP_DIAMOND_EXACT permutation=%o time=%o\n",sig,Cputime(tt);
base:=Divisor(cs[basenum]);invs:=ClassGroupAbelianInvariants(AF:ReductionDivisor:=base,Proof:=true);MDProvenClassGroupGenerators(AF);
G,_,phi,wit:=MDClassGroupDirectRRCheckedV4FastFB(AF,"C65/logs/quotient9/finite_RR.log");
hc:=[phi(Divisor(c)-Degree(c)*base):c in cs];alphac:=[hc[sig[i]]-hc[i]-Degree(cs[i])*(hc[sig[basenum]]-hc[basenum]):i in [1..#cs]];H:=sub<G|alphac>;Jc:=sub<G|hc>;ord:=&*[v:v in Invariants(G)|v ne 0];assert ord eq 21504;
widths:=[-Valuation(jAF,c):c in cs];sf:=[1..#cs];ww:=[i:i in [1..#cs]|widths[i] eq 65];assert #ww eq 2;sf[ww[1]]:=ww[2];sf[ww[2]]:=ww[1];
Ffree:=FreeAbelianGroup(#cs);beta:=hom<Ffree->H|[H!v:v in alphac]>;
rawA5:=[];
for i in [1..#cs]do v:=Ffree.(sig[sf[i]])+3*Ffree.(sf[i])-Ffree.(sig[i])-3*Ffree.i;Append(~rawA5,v);end for;
a5h:=[rawA5[i]-Degree(cs[i])*rawA5[basenum]:i in [1..#cs]];
alphaA5:=[&+[G|Eltseq(a5h[i])[j]*alphac[j]:j in [1..#cs]]:i in [1..#cs]];
for rel in Generators(Kernel(beta))do assert &+[G|Eltseq(Ffree!rel)[j]*alphaA5[j]:j in [1..#cs]] eq G!0;end for;
function ImageA5(v)
 r:=Eltseq(v@@beta);return &+[G|r[j]*alphaA5[j]:j in [1..#cs]];
end function;
aop:=hom<H->G|[ImageA5(H.i):i in [1..NumberOfGenerators(H)]]>;HK:=Kernel(aop);
Write("C65/data/quotient9/finite_group_A3_data.txt",Sprint([*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Degree(c):c in cs],widths,sig,sf,[Eltseq(v):v in hc],[Eltseq(v):v in alphac],Invariants(Jc),Invariants(H),[Eltseq(G!v):v in Generators(HK)],[Eltseq(v):v in alphaA5]*]):Overwrite:=true);
check6:=[<i,6*alphac[i] in HK>:i in[1..#cs]|widths[i] eq 65];printf "NONIC65_NONGLOBAL_SEXTIC_CUSP_TESTS %o Horder=%o Korder=%o\n",check6,#H,#HK;Write("C65/data/quotient9/non_global_cusp_tests.txt",Sprint([*basenum,check6*]):Overwrite:=true);fprintf wit,"]\n";Flush(wit);
printf "NONIC65_FINITE_ALPHA_A3_KERNEL G=%o cuspgroup=%o alpha_cusps=%o A3kernel=%o widths=%o frob5=%o time=%o\n",Invariants(G),Invariants(Jc),Invariants(H),Invariants(HK),widths,sf,Cputime(tt);print "NONIC65_FINITE_A3_COMPLETE";quit;
