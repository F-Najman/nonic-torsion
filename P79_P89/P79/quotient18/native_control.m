SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);tt:=Cputime();
load "P79_P89/P79/quotient18/finite_common.m";
load "P79_P89/data/P79/quotient18/source_places_data.m";
load "P79_P89/data/P79/quotient18/target_places1_jbase_data.m";
su:=DecAF(SigmaUData);sv:=DecAF(SigmaVData);if AFBaseCoordinate eq 1 then sig0:=hom<K->AF|su>;sig:=hom<AF->AF|sig0,sv>;else sig0:=hom<K->AF|sv>;sig:=hom<AF->AF|sig0,su>;end if;assert sig(jAF) eq jAF;assert sig(sig(sig(g))) eq g and sig(sig(sig(AF!t))) eq AF!t;
printf "NONIC79_SIGMA_NATIVE_CONTROL order3=true fixes_j=true time=%o\n",Cputime(tt);
function DecodePlace(gg)
 aa:=DecAF(gg[1]);bb:=DecAF(gg[2]);pls:=[pl:pl in Zeros(aa)|Valuation(bb,pl) gt 0];assert #pls eq 1;return pls[1];
end function;
reps:=[o[1]:o in SourceSigmaOrbits|SourceJIndices[o[1]] eq 1];k:=reps[1];kk:=[k,SourceSigmaPermutation[k],SourceSigmaPermutation[SourceSigmaPermutation[k]]];pp:=[DecodePlace(SourcePlaceGenerators[h]):h in kk];qq:=[DecodePlace(TargetPlaceGenerators[h]):h in [1,TargetSigmaPermutation[1]]];B:=qq[2]-qq[1]+pp[3]-5*pp[2]+4*pp[1];assert Degree(B) eq 0;
printf "NONIC79_RR_NATIVE_CONTROL_START fiber=1 source=%o target=1 convention=1 time=%o\n",k,Cputime(tt);dm:=Dimension(B);assert dm eq 0;printf "NONIC79_RR_NATIVE_CONTROL_COMPLETE dimension=0 time=%o\n",Cputime(tt);quit;
