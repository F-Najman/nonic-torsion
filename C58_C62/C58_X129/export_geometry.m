SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt:=Cputime();AttachSpec("mdmagma/v2/mdmagma.spec");
function Counts58(f,p)
 aa:=Reverse(Coefficients(f));ss:=[];bb:=[];
 for n in [1..9] do sn:=-n*aa[n+1]-&+[Integers()|aa[j+1]*ss[n-j]:j in [1..n-1]];Append(~ss,sn);v:=p^n+1-sn-&+[Integers()|d*bb[d]:d in Divisors(n)|d lt n];assert v mod n eq 0;Append(~bb,v div n);end for;return bb;
end function;
f29:=FrobeniusPolynomial(JOne(29),3);f58:=FrobeniusPolynomial(JOne(58),3);assert Degree(f29) eq 44 and Degree(f58) eq 156;a29:=Counts58(f29,3);a58:=Counts58(f58,3);n29:=a29;n29[1]-:=14;n58:=a58;n58[1]-:=28;assert n29 eq [0,0,14,28,14,133,280,770,2282] and &and[n58[i] eq 0:i in [1..5]];
X:=MDX1(29,GF(3));cf:=DefiningPolynomial(Curve(X));AA<x,y>:=AffineSpace(GF(3),2);R:=CoordinateRing(AA);if Rank(Parent(cf)) eq 3 then fxy:=Evaluate(cf,[R.1,R.2,R!1]);else fxy:=Evaluate(cf,[R.1,R.2]);end if;
out:=Open("C58_C62/data/C58_X129/geometry58.json","w");fprintf out,"{\"N29_frobenius_coefficients\":%o,\"N58_frobenius_coefficients\":%o,\"N29_noncuspidal_counts\":%o,\"N58_noncuspidal_counts\":%o,\"N29_affine_exponents\":%o,\"N29_affine_coefficients\":%o,\"complete\":true}\n",Coefficients(f29),Coefficients(f58),n29,n58,[Exponents(v):v in Monomials(fxy)],[Integers()!v:v in Coefficients(fxy)];delete out;
printf "NONIC58_GEOMETRY_COUNTS_COMPLETE source=%o target=%o time=%o\n",n58,n29,Cputime(tt);quit;
