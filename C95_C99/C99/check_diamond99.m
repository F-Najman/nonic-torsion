SetColumns(0);SetQuitOnError(true);
a:=2;N:=99;o:=1;while a^o mod N notin {1,N-1} do o+:=1;end while;assert o eq 15;
print "NONIC_C99_DIAMOND_PASS a=2 projective_order=15 complete=true";quit;
