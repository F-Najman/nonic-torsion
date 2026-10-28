SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
load "P71/quotient26/finite_common.m";
load "P71/data/quotient26/source_places_data.m";
su:=DecAF(SigmaUData);sv:=DecAF(SigmaVData);if AFBaseCoordinate eq 1 then sig0:=hom<K->AF|su>;sig:=hom<AF->AF|sig0,sv>;else sig0:=hom<K->AF|sv>;sig:=hom<AF->AF|sig0,su>;end if;
assert sig(jAF) eq jAF;uu:=uAF;vv:=vAF;for _ in[1..5]do uu:=sig(uu);vv:=sig(vv);end for;assert uu eq uAF and vv eq vAF;
printf "NONIC71_SIGMA_NATIVE_CONTROL order5=true fixes_j=true time=%o\n",Cputime(tt);
function DecodePlace(gg)
 aa:=DecAF(gg[1]);bb:=DecAF(gg[2]);pls:=[pl:pl in Zeros(aa)|Valuation(bb,pl) gt 0];assert #pls eq 1;return pls[1];
end function;
load "P71/data/quotient26/target_places1_jbase_data.m";
k:=SourceSigmaOrbits[1][1];pp:=[DecodePlace(SourcePlaceGenerators[i]):i in[k,SourceSigmaPermutation[k]]];q:=DecodePlace(TargetPlaceGenerators[1]);B:=q-pp[2]-3*pp[1];assert Degree(B) eq 0;
printf "NONIC71_RR_NATIVE_CONTROL_START source=%o target_fiber1 convention1 time=%o\n",k,Cputime(tt);assert Dimension(B) eq 0;printf "NONIC71_RR_NATIVE_CONTROL_COMPLETE degree7 convention1 dimension0 time=%o\n",Cputime(tt);
load "P71/data/quotient26/target_places2_jbase_data.m";
q1:=DecodePlace(TargetPlaceGenerators[1]);
load "P71/data/quotient26/target_places3_jbase_data.m";
q2:=DecodePlace(TargetPlaceGenerators[1]);k:=SourceSigmaOrbits[2][1];pp:=[DecodePlace(SourcePlaceGenerators[i]):i in[k,SourceSigmaPermutation[k]]];B:=q1+q2-pp[1]-3*pp[2];assert Degree(B) eq 0;
printf "NONIC71_RR_NATIVE_CONTROL_START source=%o target_fibers2and3 convention2 time=%o\n",k,Cputime(tt);assert Dimension(B) eq 0;printf "NONIC71_RR_NATIVE_CONTROL_COMPLETE degree8 convention2 dimension0 time=%o\n",Cputime(tt);print "NONIC71_NATIVE_CONTROLS_COMPLETE";quit;
