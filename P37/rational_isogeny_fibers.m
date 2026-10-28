SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
for lab in ["1225h1","1225h2"] do
 E:=EllipticCurve(lab);assert jInvariant(E) in [-9317,-162677523113838677];
 fs:=Factorization(DivisionPolynomial(E,37));ds:=[<Degree(f[1]),f[2]>:f in fs];
 assert &and[IsEven(z[1]) and z[2] eq 1:z in ds];
 printf "NONIC37_RATIONAL_ISOGENY_FIBER label=%o a=%o j=%o x_factor_degrees=%o all_even=true time=%o\n",lab,aInvariants(E),jInvariant(E),ds,Cputime(tt);
end for;
print "NONIC37_RATIONAL_ISOGENY_FIBERS_COMPLETE";quit;
