SetColumns(0);SetQuitOnError(true);N:=95;a:=11;
assert a^3 mod N eq 1 and a ne 1;
chars:=[chi:chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N))|IsEven(chi) and Order(chi) eq 12];
assert #chars eq 1 and &and[Evaluate(chi,a) eq 1:chi in chars];
print "NONIC_95_FINITE_DIAMOND_IMAGE_KERNEL_PASS a=11 order=3 bad_character_order=12 complete=true";quit;
