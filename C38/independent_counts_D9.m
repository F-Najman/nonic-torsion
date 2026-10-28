SetColumns(0);SetMemoryLimit(2*10^9);t0:=Cputime();
AttachSpec("gl2/magma.spec");
for tup in [<35,3>,<38,3>,<44,3>,<46,3>,<50,3>] do
 n,p:=Explode(tup);GG:=GL(2,Integers(n));
 HH:=sub<GG|[GG![1,0,1,1]] cat [GG![a,0,0,1]:a in [1..n]|GCD(a,n) eq 1 ]>;
 zs:=GL2PointCounts(HH,p,9);
 ys:=[zs[d]-GL2RationalCuspCount(HH,p^d):d in [1..9]];
 assert &+GL2jCounts(HH,p) eq ys[1];closed:=[];
 for d in [1..9] do
  z:=ys[d]-&+([0] cat [e*closed[e]:e in [1..d-1]|d mod e eq 0]);
  assert z mod d eq 0;Append(~closed,z div d);
 end for;
 printf "DEGREE9_FAST_COUNT m=1 n=%o p=%o genus=%o rational_point_counts=%o closed_place_counts=%o time=%o\n",n,p,GL2Genus(HH),ys,closed,Cputime(t0);
end for;
print "DEGREE9_FAST_COUNTS_DONE complete=true",Cputime(t0);
quit;
