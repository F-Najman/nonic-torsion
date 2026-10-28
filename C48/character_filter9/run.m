SetColumns(0);SetQuitOnError(true);SetMemoryLimit(6*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");load "C48/data/portable_inputs.m";load "C48/data/character_filter9/functions.m";load "C48/data/parameter9/records.m";assert Complete and #Records eq 13513;
Y:=MDX11(2,24,GF(5));FC:=FunctionField(Curve(Y));F0,am:=AlgorithmicFunctionField(FC);F:=RationalExtensionRepresentation(F0);B:=BaseField(F);P:=PolynomialRing(GF(5));assert [F!am(FC.i):i in [1..2]] eq [F.1,F!B.1];
function Rat(v)return B!(P!v[1])/B!(P!v[2]);end function;
assert Coefficients(DefiningPolynomial(F)) eq [Rat(v):v in PortableFieldEquation];
for data in FilterFunctions do f:=F![Rat(v):v in data];D:=Divisor(f);ps,es:=Support(D);assert &and[e mod 2 eq 0:e in es];end for;
CM:=Matrix(GF(2),[[Integers()!x:x in v]:v in FilterCharacterRows]);HM:=Matrix(GF(2),[[Integers()!x:x in v[1..8]]:v in PortableRationalCusps]);KM:=Matrix(GF(2),[[Integers()!x:x in v]:v in FilterCombinations]);assert Rank(CM) eq 8 and Rank(KM*CM) eq 4 and IsZero(KM*CM*Transpose(HM));
R:=PolynomialRing(GF(5));K<z>:=ext<GF(5)|R!FieldPolynomial>;
function UV(b,c)
 if 1-c eq 0 then return false,[];end if;t:=2*b/(1-c)-1;if t eq 0 then return false,[];end if;q:=((1-c)*t^2-2)/(2*t);den:=t*(q^2+1)+2*q;if den*q eq 0 then return false,[];end if;ss:=t*(q^2-1)/den;dd:=(ss+2)/q;return true,[(ss+dd)/2,(ss-dd)/2];
end function;
needed:=[i:i in [1..8]|&or[row[i] ne 0:row in FilterCombinations]];assert needed eq [2,3,6,7,8];
function SignAt(i,u,v)
 value:=K!0;power:=K!1;
 for pair in FilterFunctions[i] do den:=Evaluate(R!pair[2],v);if den eq 0 then return 0;end if;value+:=Evaluate(R!pair[1],v)/den*power;power*:=u;end for;
 if value eq 0 then return 0;end if;sg:=(value/K!FilterBaseValues[i])^((#K-1) div 2);assert sg in {K!1,-K!1};return sg eq 1 select 1 else -1;
end function;
out:=Open("C48/data/character_filter9/results.py","w");fprintf out,"rows=[\n";remaining:=[];excluded:=[];noSource:=0;fallback:=0;weightedRemaining:=0;
for i in [1..#Records] do
 row:=Records[i];if row[5] eq 0 then noSource+:=1;fprintf out,"[%o,-1,[]],\n",i;continue;end if;
 b:=Evaluate(R!row[2],z);c:=Evaluate(R!row[3],z);ok,uv:=UV(b,c);signs:=[0:j in [1..8]];pass:=true;
 if ok then
  u,v:=Explode(uv);assert &+[Evaluate(R!(PortableFieldEquation[j][1]),v)/Evaluate(R!(PortableFieldEquation[j][2]),v)*u^(j-1):j in [1..#PortableFieldEquation]] eq 0;
  for k in needed do signs[k]:=SignAt(k,u,v);end for;
  if &or[signs[k] eq 0:k in needed] then fallback+:=1;
  else for combo in FilterCombinations do if &*[K!(signs[j]^combo[j]):j in [1..8]|combo[j] ne 0] ne 1 then pass:=false;break;end if;end for;end if;
 else fallback+:=1;end if;
 if pass then Append(~remaining,i);weightedRemaining+:=row[5];else Append(~excluded,i);end if;
 fprintf out,"[%o,%o,%o],\n",i,pass select 1 else 0,signs;Flush(out);
 if i mod 500 eq 0 then printf "NONIC48_CHARACTER9_PROGRESS tested=%o excluded=%o remaining=%o fallback=%o time=%o\n",i,#excluded,#remaining,fallback,Cputime(tt);end if;
end for;
fprintf out,"]\nremaining_indices=%o\nexcluded_indices=%o\nno_source_orbits=%o\nfallback=%o\ncomplete=True\n",remaining,excluded,noSource,fallback;delete out;
Write("C48/data/character_filter9/remaining.m",Sprintf("RemainingIndices:=%o;",remaining):Overwrite:=true);
printf "NONIC48_CHARACTER9_DONE tested=13513 no_source=%o excluded=%o remaining=%o weighted_source_remaining=%o fallback=%o complete=true time=%o\n",noSource,#excluded,#remaining,weightedRemaining,fallback,Cputime(tt);quit;
