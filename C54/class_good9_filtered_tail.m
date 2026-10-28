/* Load after the fresh proven setup_class_queries.m. */
SetColumns(0);SetQuitOnError(true);tquery:=Cputime();
assert assigned Shard and Shard in [0..3];NShards:=4;
assert Invariants(G5) eq [3,3,342175446,0] and G5AF cmpeq G5 and #candidate2 eq 2412;
assert phi(afinv(G5.3)) eq G5.3;
load "common/pair_af_place_linear.m";
load "C54/data/quadratic_indices.m";
assert #Skip54 eq 8087 and #Remaining54 eq 8048 and #Fallback54 eq 0 and Seqset(Skip54) meet Seqset(Remaining54) eq {} and Seqset(Skip54 cat Remaining54) eq Seqset([1..16135]);
data:=eval Read("C54/data/good9_pointdata.txt");def,records:=Explode(data);assert #records eq 16135;
rr:=PolynomialRing(GF(5));ff:=rr!def;assert IsIrreducible(ff) and Degree(ff) eq 9;K9<a9>:=ext<GF(5)|ff>;
skips:=[i:i in Skip54|(i-1) mod 4 eq Shard];todo:=[i:i in Remaining54|(i-1) mod 4 eq Shard];
if not assigned OutFile then OutFile:=Sprintf("C54/data/filtered9_classes_shard%o.py",Shard);end if;
out:=Open(OutFile,"w");fprintf out,"shard=%o\nnshards=4\nquadratic_skip_indices=%o\nrelations2=%o\nrelations5=%o\ncusp2=%o\ncusp5=%o\ncusp_signatures=%o\ncandidate2=%o\ntests=[\n",Shard,skips,[Eltseq(r):r in Rows(RelationMatrix(G2))],[Eltseq(r):r in Rows(RelationMatrix(G5))],[Eltseq(v):v in qc2],[Eltseq(v):v in qc5],sig2,[Eltseq(v):v in candidate2];Flush(out);
tested:=0;members:=0;passing:=0;
for i in todo do
 xyco,n,r2,jco,nr3,orbitLength:=Explode(records[i]);assert n mod 54 eq 0 and r2 in [1,3] and nr3 eq 2 and orbitLength in [1,3,9];
 xy:=[Evaluate(rr!c,a9):c in xyco];assert LCM([Degree(MinimalPolynomial(z,GF(5))):z in xy]) eq 9;
 Q:=C48PairAFPlaceLinear(coords,xy);assert Degree(Q) eq 9;aa:=phi(Divisor(Q)-9*A0);member:=aa in H5;passes:=false;im2:=[];preimage:=[];
 if member then
  members+:=1;v:=aa@@h5;assert h5(v) eq aa;im2:=Eltseq(h2(v));preimage:=Eltseq(v);if h2(v) in candidate2 then passes:=true;passing+:=1;end if;
 end if;
 tested+:=1;fprintf out,"[%o,%o,%o,%o,%o,%o],\n",i,Eltseq(aa),member select 1 else 0,passes select 1 else 0,im2,preimage;Flush(out);
 if tested le 3 or tested mod 25 eq 0 then printf "NONIC54_FILTERED9_CLASS_PROGRESS shard=%o tested=%o members=%o passing=%o query_cpu=%o\n",Shard,tested,members,passing,Cputime(tquery);end if;
 if assigned MaxTests and tested ge MaxTests then fprintf out,"]\nbenchmark_only=True\n";delete out;printf "NONIC54_FILTERED9_BENCHMARK shard=%o tested=%o query_cpu=%o complete=false\n",Shard,tested,Cputime(tquery);quit;end if;
end for;
fprintf out,"]\ncomplete=True\n";delete out;assert tested eq #todo;
printf "NONIC54_FILTERED9_CLASS_DONE shard=%o class_queries=%o quadratic_skips=%o members=%o passing=%o total_records=%o complete=true query_cpu=%o\n",Shard,tested,#skips,members,passing,tested+#skips,Cputime(tquery);quit;
