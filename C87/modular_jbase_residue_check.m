SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt:=Cputime();
load "C87/data/finite_j_model_snapshot.m";
load "C87/data/modular_jbase_residue_data.m";
load "C87/data/jbase_relation_modular_data.m";
PF<t>:=PolynomialRing(GF(2));assert #ModularBasePrimes eq 81 and #Seqset(ModularBasePrimes)eq 81 and #ModularNormResidues eq 81;assert JBaseDegreeJ eq 31 and JBaseDegreeParameter eq 480;
function CoefficientValue(c,a)
 n:=Evaluate(PF!c[1],a);d:=Evaluate(PF!c[2],a);assert d ne 0;return n/d;
end function;
function BinaryValue(v,a)
 r:=Parent(a)!0;b:=Parent(a)!1;
 for _ in[0..11]do if IsOdd(v)then r+:=b;end if;v:=v div 2;b*:=a;end for;assert v eq 0;return r;
end function;
for idx in[1..81]do
 v:=ModularBasePrimes[idx];mp:=PF![(v div 2^k)mod 2:k in[0..12]];assert Degree(mp)eq 12 and IsIrreducible(mp);F<a>:=ext<GF(2)|mp>;
 ff:=[CoefficientValue(c,a):c in AFPolynomialData];jj:=[CoefficientValue(c,a):c in JFunctionData];assert #ff eq 32 and #jj eq 31 and ff[32]ne 0;top:=ff[32];ff:=[b/top:b in ff];cols:=[];vv:=jj;
 for k in[1..31]do Append(~cols,vv);top:=vv[31];vv:=[F!0]cat vv[1..30];vv:=[vv[i]+top*ff[i]:i in[1..31]];end for;
 A:=Transpose(Matrix(F,cols));cp:=CharacteristicPolynomial(A);assert Degree(cp)eq 31 and Coefficient(cp,30)eq Trace(A)and Coefficient(cp,0)eq Determinant(A);
 expected:=[BinaryValue(c,a):c in ModularNormResidues[idx]];assert Eltseq(cp)eq expected;
 relation:=[Evaluate(PF!c,a):c in JBasePolynomialCoefficients];assert relation[32]ne 0;top:=relation[32];assert [c/top:c in relation]eq expected;
 if idx mod 10 eq 0 then printf "NONIC87_NATIVE_MODULAR_NORM_PROGRESS primes=%o time=%o\n",idx,Cputime(tt);end if;
end for;
print "NONIC87_NATIVE_MODULAR_NORM_CHECK_COMPLETE primes81 modulusdegree972 heightbound480 all_exact_residues_agree=true";quit;
