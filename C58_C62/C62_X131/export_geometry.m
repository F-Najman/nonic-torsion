SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt:=Cputime();AttachSpec("mdmagma/v2/mdmagma.spec");
function Counts62(f,p)
 aa:=Reverse(Coefficients(f));ss:=[];bb:=[];
 for n in [1..9] do sn:=-n*aa[n+1]-&+[Integers()|aa[j+1]*ss[n-j]:j in [1..n-1]];Append(~ss,sn);v:=p^n+1-sn-&+[Integers()|d*bb[d]:d in Divisors(n)|d lt n];assert v mod n eq 0;Append(~bb,v div n);end for;return bb;
end function;
f31:=FrobeniusPolynomial(JOne(31),3);f62:=FrobeniusPolynomial(JOne(62),3);assert Degree(f31) eq 52 and Degree(f62) eq 182;a31:=Counts62(f31,3);a62:=Counts62(f62,3);n31:=a31;n31[1]-:=15;n62:=a62;n62[1]-:=30;assert &and[n62[i] eq 0:i in [1..4]];
X:=MDX1(31,GF(3));cf:=DefiningPolynomial(Curve(X));AA<x,y>:=AffineSpace(GF(3),2);R:=CoordinateRing(AA);if Rank(Parent(cf)) eq 3 then fxy:=Evaluate(cf,[R.1,R.2,R!1]);else fxy:=Evaluate(cf,[R.1,R.2]);end if;
out:=Open("C58_C62/data/C62_X131/geometry62.json","w");fprintf out,"{\"N31_frobenius_coefficients\":%o,\"N62_frobenius_coefficients\":%o,\"N31_noncuspidal_counts\":%o,\"N62_noncuspidal_counts\":%o,\"N31_affine_exponents\":%o,\"N31_affine_coefficients\":%o,\"complete\":true}\n",Coefficients(f31),Coefficients(f62),n31,n62,[Exponents(v):v in Monomials(fxy)],[Integers()!v:v in Coefficients(fxy)];delete out;
printf "NONIC62_GEOMETRY_COUNTS_COMPLETE source=%o target=%o time=%o\n",n62,n31,Cputime(tt);quit;
