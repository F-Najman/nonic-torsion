/* Portable candidate presentation, saved BEFORE any divisor-class queries.
   Proof=false is only a relation-search mode. The output is NOT a proven
   class group until every presentation relation and prime-order line passes
   the independent RR verifier and the Frobenius order agrees. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
assert assigned N and assigned p and assigned OutDir;
if Type(N) eq MonStgElt then N:=StringToInteger(N);end if;
if Type(p) eq MonStgElt then p:=StringToInteger(p);end if;
assert IsPrime(p) and N mod p ne 0;
if not assigned LogDir then LogDir:=OutDir;end if;log:=Open(LogDir cat "/build.log","w");
procedure Note(s) fprintf log,"%o\n",s;Flush(log);print s;end procedure;
X:=MDX1(N,GF(p));C:=Curve(X);FF:=FunctionField(C);AF,am:=AlgorithmicFunctionField(FF);
AF:=RationalExtensionRepresentation(AF);xy0:=[AF!am(FF.i):i in [1..2]];
K:=BaseField(AF);assert Type(K) eq FldFunRat and #ConstantField(K) eq p;
R<z>:=PolynomialRing(K);f0:=R!DefiningPolynomial(AF);f0/:=LeadingCoefficient(f0);d:=Degree(f0);
den:=LCM([Denominator(a):a in Coefficients(f0)]);
f:=&+[Coefficient(f0,i)*(K!den)^(d-i)*z^i:i in [0..d]];
assert LeadingCoefficient(f) eq 1 and &and[Degree(Denominator(a)) eq 0:a in Coefficients(f)];
F<u>:=FunctionField(f);
fw:=hom<AF->F|hom<K->F|F!K.1>,u/(F!den)>;
bw:=hom<F->AF|hom<K->AF|AF!K.1>,(AF!den)*AF.1>;
assert fw(bw(F.1)) eq F.1 and bw(fw(AF.1)) eq AF.1;
coords:=[fw(a):a in xy0];assert RationalExtensionRepresentation(F) cmpeq F;
g:=Genus(F);assert g eq Genus(X);
frob:=FrobeniusPolynomial(JOne(N),p);known:=Integers()!Evaluate(frob,1);
if assigned ExpectedOrder then assert known eq ExpectedOrder;end if;
Note(Sprintf("PRESENTATION_MODEL N=%o p=%o genus=%o integral_monic_degree=%o natural_pole_degree=%o known_Pic_order=%o time=%o",N,p,g,d,Degree(Denominator(Divisor(F!K.1))),known,Cputime(tt)));
inv:=ClassGroupAbelianInvariants(F:Proof:=false);
assert assigned F`cld`rels and assigned F`cld`fb;
fb:=[a`place:a in F`cld`fb];rels:=[Eltseq(v):v in F`cld`rels];
Note(Sprintf("PRESENTATION_RELATIONS_READY invariants=%o fb_size=%o relation_rows=%o fb_bound=%o native_proven=%o time=%o",inv,#fb,#rels,F`cld`fbbound,F`cld`proven,Cputime(tt)));
function EncRat(a)
 return [[Integers()!c:c in Coefficients(Numerator(a))],[Integers()!c:c in Coefficients(Denominator(a))]];
end function;
function EncElt(a)
 return [EncRat(c):c in Eltseq(F!a)];
end function;
function DivVector(D,fb)
 ps,es:=Support(D);v:=[0:i in [1..#fb]];
 for i in [1..#ps] do j:=Index(fb,ps[i]);assert j gt 0;v[j]+:=es[i];end for;
 return v;
end function;
A1:=F`cld`A1;assert Degree(A1) eq 1;Ainf:=Denominator(Divisor(F!K.1));
out:=Open(OutDir cat "/presentation.json","w");
fprintf out,"{\"status\":\"candidate_presentation_only\",\"N\":%o,\"p\":%o,\"genus\":%o,\"field_degree\":%o,\"field_polynomial\":%o,\"modular_coordinates\":%o,\"known_Pic_order\":%o,\"frobenius_coefficients\":%o,\"candidate_invariants\":%o,\"native_proven\":%o,\"factor_basis_bound\":%o,\"relations\":%o,\"degree_one_vector\":%o,\"natural_pole_vector\":%o,\"factor_basis\":[",N,p,g,d,[EncRat(a):a in Coefficients(DefiningPolynomial(F))],[EncElt(a):a in coords],known,Coefficients(frob),inv,F`cld`proven,F`cld`fbbound,rels,DivVector(A1,fb),DivVector(Ainf,fb);
for i in [1..#fb] do
 P:=fb[i];I:=Ideal(P);bb:=[F!a:a in Basis(I)];
 if i gt 1 then fprintf out,",";end if;
 fprintf out,"{\"index\":%o,\"degree\":%o,\"finite\":%o,\"ideal_basis\":%o}",i,Degree(P),IsFinite(P),[EncElt(a):a in bb];Flush(out);
end for;
fprintf out,"],\"checkpoint_complete\":true}\n";delete out;
Note(Sprintf("PRESENTATION_CHECKPOINT_SAVED N=%o fb_size=%o relation_rows=%o no_place_mapping_started=true time=%o",N,#fb,#rels,Cputime(tt)));
delete log;quit;
