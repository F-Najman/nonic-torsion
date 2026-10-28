SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);
AttachSpec("gl2/magma.spec");
n:=63;p:=2;GG:=GL(2,Integers(n));
HH:=sub<GG|[GG![1,0,1,1]] cat [GG![a,0,0,1]:a in [1..n]|GCD(a,n) eq 1]>;
zs:=GL2PointCounts(HH,p,9);ys:=[zs[d]-GL2RationalCuspCount(HH,p^d):d in [1..9]];closed:=[];
for d in [1..9] do
 v:=ys[d]-&+([0] cat [e*closed[e]:e in [1..d-1]|d mod e eq 0]);assert v mod d eq 0;Append(~closed,v div d);
end for;
assert GL2Genus(HH) eq 97 and closed eq [0,0,0,0,0,0,54,54,72];
printf "NONIC63_SOURCE_COUNTS genus=97 rational=%o closed=%o complete=true\n",ys,closed;
quit;
