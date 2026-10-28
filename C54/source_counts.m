SetColumns(0);SetQuitOnError(true);tt:=Cputime();
AttachSpec("gl2/magma.spec");
N:=54;p:=5;GG:=GL(2,Integers(N));HH:=sub<GG|[GG![1,0,1,1]] cat [GG![a,0,0,1]:a in [1..N]|GCD(a,N) eq 1]>;
z:=GL2PointCounts(HH,p,9);y:=[z[d]-GL2RationalCuspCount(HH,p^d):d in [1..9]];closed:=[];
for d in [1..9] do a:=y[d]-&+([0] cat [e*closed[e]:e in [1..d-1]|d mod e eq 0]);assert a mod d eq 0;Append(~closed,a div d);end for;
assert GL2Genus(HH) eq 52 and closed[1..8] eq [0,0,54,261,612,3033,11412,48087];
Write("C54/data/source_counts_data.m",Sprintf("Source54Counts:=%o;\n",closed):Overwrite:=true);
printf "NONIC54_SOURCE_COUNTS genus=52 p=5 rational=%o closed=%o complete=true time=%o\n",y,closed,Cputime(tt);quit;
