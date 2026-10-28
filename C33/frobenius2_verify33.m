SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt:=Cputime();
load "C33/data/frobenius2_input33.m";
f:=FrobeniusPolynomial(JOne(33),2);assert Coefficients(f)eq ExpectedCoefficients and Degree(f)eq 42 and Evaluate(f,1)eq 635604750;
printf "NONIC33_FRESH_FROBENIUS2_VERIFIED degree42_order635604750=true all_coefficients_exact=true time=%o\n",Cputime(tt);quit;
