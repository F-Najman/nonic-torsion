SetColumns(0);SetQuitOnError(true);
N:=841;a:=12;o:=1;while a^o mod N notin {1,N-1} do o+:=1;end while;assert o eq 58;
chars:=[chi:chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N))|IsEven(chi) and Order(chi) in {1,7}];
assert &and[Evaluate(chi,a) eq 1:chi in chars];
print "NONIC_841_DIAMOND_KERNEL_PASS order=58 bad_character_orders=1,7 exact=true";quit;
