SetColumns(0);SetQuitOnError(true);SetMemoryLimit(32*10^9);tt:=Cputime();
load "P79_P89/P89/quotient27/finite_common.m";
load "P79_P89/data/P89/quotient27/source_places_data.m";
assert assigned FiberIndex;i:=StringToInteger(FiberIndex);assert i in [1..3];
load "P79_P89/data/P89/quotient27/target_places1_jv_data.m";
function DecodePlace(gg)
 aa:=DecAF(gg[1]);bb:=DecAF(gg[2]);pls:=[pl:pl in Zeros(aa)|Valuation(bb,pl) gt 0];assert #pls eq 1;return pls[1];
end function;

k:=5;kk:=[k];for _ in [1..3]do Append(~kk,SourceSigmaPermutation[kk[#kk]]);end for;
pp:=[DecodePlace(SourcePlaceGenerators[h]):h in kk];qq:=[DecodePlace(TargetPlaceGenerators[h]):h in [1,TargetSigmaPermutation[1]]];
B:=qq[2]-qq[1]-pp[4]+pp[3]-5*pp[2]+5*pp[1];assert Degree(B) eq 0;
printf "NONIC89_RR_DIMENSION_PILOT_START degree=0 time=%o\n",Cputime(tt);
dim:=Dimension(B);printf "NONIC89_RR_DIMENSION_PILOT_RESULT dimension=%o time=%o\n",dim,Cputime(tt);quit;
