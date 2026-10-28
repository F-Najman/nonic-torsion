SetColumns(0);SetQuitOnError(true);SetMemoryLimit(32*10^9);tt:=Cputime();
load "P79_P89/P79/quotient18/finite_common.m";
pls:=Places(AF,9);eligible:=[];jindices:=[];zs:=[* *];keys:=[* *];
for pl in pls do
 if Valuation(jAF,pl) lt 0 then continue;end if;
 jv:=Evaluate(jAF,pl);jj:=[Integers()!c:c in Eltseq(MinimalPolynomial(jv,GF(2)))];
 ix:=Index([s:s in SourceJPolynomials],jj);if ix eq 0 then continue;end if;
 Append(~eligible,pl);Append(~jindices,ix);z:=PlaceCoordinates(pl,Z);Append(~zs,z);Append(~keys,JKey(z,jv));
end for;
assert #eligible eq #Seqset([k:k in keys]);perm:=[];
for k in [1..#eligible] do
 jv:=Evaluate(jAF,eligible[k]);key:=JKey(SigmaCoordinates(zs[k]),jv);ii:=Index(keys,key);assert ii ne 0;Append(~perm,ii);
end for;
assert Sort(perm) eq [1..#perm];for i in [1..#perm]do assert perm[perm[perm[i]]] eq i;end for;
seen:={};orbs:=[];for i in [1..#eligible]do if i in seen then continue;end if;orb:=[i];k:=perm[i];while k ne i do Append(~orb,k);k:=perm[k];end while;Append(~orbs,orb);seen join:=Seqset(orb);end for;
gens:=[];for pl in eligible do a,b:=TwoGenerators(pl);Append(~gens,[CodeAF(a),CodeAF(b)]);end for;
out:="P79_P89/data/P79/quotient18/source_places_data.m";
PrintFile(out,Sprintf("SourcePlaceGenerators:=%o;SourcePlaceKeys:=%o;SourceJIndices:=%o;SourceSigmaPermutation:=%o;SourceSigmaOrbits:=%o;",gens,keys,jindices,perm,orbs):Overwrite:=true);
PrintFile(out,"print \"NONIC79_QUOTIENT_SOURCE_PLACES_READY\";");
printf "NONIC79_QUOTIENT_SOURCE_PLACES all_degree9=%o eligible=%o sigma_orbits=%o j_multiplicities=%o time=%o\n",#pls,#eligible,[#o:o in orbs],[#[i:i in jindices|i eq k]:k in [1..2]],Cputime(tt);quit;
