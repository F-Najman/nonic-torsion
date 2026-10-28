SetColumns(0);SetQuitOnError(true);
N:=899;a:=3;o:=1;while a^o mod N notin {1,N-1} do o+:=1;end while;assert o eq 420;
chars:=[c:c in GaloisConjugacyRepresentatives(FullDirichletGroup(N))|IsEven(c) and Order(c) in {1}];
assert &and[Evaluate(c,a) eq 1:c in chars];
print "NONIC_899_DIAMOND_KERNEL_PASS order=420 exact=true";quit;
