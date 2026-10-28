SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);tt:=Cputime();
load "P73/quotient19/finite_common.m";
load "P73/data/quotient19/source_places_data.m";
load "P73/data/quotient19/target_places1_jbase_data.m";
su:=DecAF(SigmaUData);sv:=DecAF(SigmaVData);if AFBaseCoordinate eq 1 then sig0:=hom<K->AF|su>;sig:=hom<AF->AF|sig0,sv>;else sig0:=hom<K->AF|sv>;sig:=hom<AF->AF|sig0,su>;end if;assert sig(jAF) eq jAF;assert sig(sig(sig(sig(g)))) eq g and sig(sig(sig(sig(AF!t)))) eq AF!t;
printf "NONIC73_SIGMA_NATIVE_CONTROL order4=true fixes_j=true time=%o\n",Cputime(tt);
function DecodePlace(gg)
 aa:=DecAF(gg[1]);bb:=DecAF(gg[2]);pls:=[pl:pl in Zeros(aa)|Valuation(bb,pl) gt 0];assert #pls eq 1;return pls[1];
end function;
k:=SourceSigmaOrbits[1][1];kk:=[k];for _ in [1..3]do Append(~kk,SourceSigmaPermutation[kk[#kk]]);end for;pp:=[DecodePlace(SourcePlaceGenerators[h]):h in kk];qq:=[DecodePlace(TargetPlaceGenerators[h]):h in [1,TargetSigmaPermutation[1]]];B:=qq[2]-qq[1]-pp[4]+pp[3]-3*pp[2]+3*pp[1];assert Degree(B) eq 0;
printf "NONIC73_RR_NATIVE_CONTROL_START source=%o target=1 convention=1 time=%o\n",k,Cputime(tt);dm:=Dimension(B);assert dm eq 0;printf "NONIC73_RR_NATIVE_CONTROL_COMPLETE dimension=0 time=%o\n",Cputime(tt);quit;
