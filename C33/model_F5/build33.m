SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
X:=MDX1(33,GF(5));C:=Curve(X);assert Genus(X)eq 21;os:=CuspOrbitsQ(X);cs:=Cusps(X);
assert #os eq 17 and Sort([Degree(o):o in os])eq [1:i in[1..10]]cat[2:i in[1..5]]cat[10,10];
assert #Seqset([CuspSignature(X,Support(o)[1]):o in os])eq 17;
F:=AlgorithmicFunctionField(FunctionField(C));assert Genus(F)eq 21;K:=BaseField(F);k:=ConstantField(F);
function CodeRat(a)return [[Integers()!c:c in Eltseq(Numerator(a))],[Integers()!c:c in Eltseq(Denominator(a))]];end function;
function CodeF(a)return [CodeRat(c):c in Eltseq(F!a)];end function;
out:=Open("C33/data/model_F5/model_cusps33.json","w");
fprintf out,"{\"N\":33,\"p\":5,\"genus\":21,\"AF_degree\":%o,\"field_polynomial\":%o,\"global_cusp_blocks\":[",Degree(F),[CodeRat(c):c in Coefficients(DefiningPolynomial(F))];
for i in[1..#os]do
 o:=os[i];sig:=CuspSignature(X,Support(o)[1]);assert &and[CuspSignature(X,P)eq sig:P in Support(o)];DD:=FunctionFieldDivisor(o);ps,es:=Support(DD);
 assert Degree(DD)eq Degree(o);
 if i gt 1 then fprintf out,",";end if;
 fprintf out,"{\"signature\":%o,\"global_degree\":%o,\"places\":[",sig,Degree(o);
 for j in[1..#ps]do
  if j gt 1 then fprintf out,",";end if;P:=ps[j];
  fprintf out,"{\"degree\":%o,\"multiplicity\":%o,\"finite\":%o,\"ideal_basis\":%o}",Degree(P),es[j],IsFinite(P)select 1 else 0,[CodeF(v):v in Basis(Ideal(P))];
 end for;fprintf out,"]}";Flush(out);
end for;fprintf out,"],\"complete\":true}\n";delete out;
printf "NONIC33_F5_MODEL_CUSPS_COMPLETE genus21=true global_blocks17=true all_signatures_unique=true AFdegree=%o time=%o\n",Degree(F),Cputime(tt);quit;
