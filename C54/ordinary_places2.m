SetColumns(0);SetQuitOnError(true);tt:=Cputime();AttachSpec("mdmagma/v2/mdmagma.spec");X:=MDX1(27,GF(2));
for d in [6..9] do ps:=NoncuspidalPlaces(X,d);ss:=[P:P in ps|jInvariant(EllipticCurve(X,P)) eq 0];printf "NONIC54_ORDINARY_COUNTS degree=%o all=%o supersingular=%o ordinary=%o time=%o\n",d,#ps,#ss,#ps-#ss,Cputime(tt);end for;
print "NONIC54_ORDINARY_COUNTS_DONE complete=true";quit;
