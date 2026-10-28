SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
load "C54/preflight_ordinary_module9.m";
FF5:=FunctionField(C5);AF5,am5:=AlgorithmicFunctionField(FF5);AF5:=RationalExtensionRepresentation(AF5);coords:=[AF5!am5(FF5.i):i in [1..2]];
G5AF,afinv,phi:=ClassGroup(AF5:Proof:=true);assert G5AF cmpeq G5;A0:=FunctionFieldDivisor(Divisor(c05));assert Degree(A0) eq 1;
QuadraticFile:="C54/data/quadratic_prepass_function.json";
load "C54/quadratic_character.m";
save "C54/data/class_preflight.workspace";
printf "NONIC54_CLASS_PREFLIGHT_NATIVE_WORKSPACE_SAVED maps_preserved=true before_all_point_queries=true time=%o\n",Cputime(tt);
data:=eval Read("C54/data/good9_pointdata.txt");def,records:=Explode(data);assert #records eq 16135;
Rp:=PolynomialRing(GF(5));pol:=Rp!def;assert IsIrreducible(pol) and Degree(pol) eq 9;KK<a9>:=ext<GF(5)|pol>;
skip:=[];remaining:=[];fallback:=[];signs:=[];
for i in [1..#records] do
 xyco,n,r2,jco,nr3,orbitLength:=Explode(records[i]);xy:=[Evaluate(Rp!c,a9):c in xyco];sign:=QuadraticCharacterXY(xy);Append(~signs,sign);
 if sign eq -1 then Append(~skip,i);else Append(~remaining,i);if sign eq 0 then Append(~fallback,i);end if;end if;
end for;
assert #skip+#remaining eq 16135;
out:=Open("C54/data/quadratic_indices.m","w");fprintf out,"Skip54:=%o;\nRemaining54:=%o;\nFallback54:=%o;\n",skip,remaining,fallback;delete out;
out:=Open("C54/data/quadratic_prepass_results.json","w");fprintf out,"{\"record_count\":16135,\"field_polynomial\":%o,\"signs\":%o,\"nonsquare_indices\":%o,\"remaining_indices\":%o,\"fallback_indices\":%o,\"complete\":true}\n",def,signs,skip,remaining,fallback;delete out;
printf "NONIC54_QUADRATIC_PREPASS_DONE records=16135 nonsquare=%o remaining=%o fallback=%o all_exact=true complete=true time=%o\n",#skip,#remaining,#fallback,Cputime(tt);quit;
