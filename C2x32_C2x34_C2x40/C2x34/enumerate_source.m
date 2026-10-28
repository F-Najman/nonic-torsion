SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");Attach("common/fast_diamond_hash_m2.m");Attach("C2x32_C2x34_C2x40/C2x34/deps/fast_full_aut_hash_m2_34.m");
X:=MDX11(2,34,GF(3));expected:=[0,0,0,36,48,48,144,738,2160];source:=[* *];small:=[* *];coverage:=[* *];
for d in [4..9] do
 xs:=NoncuspidalPlaces(X,d);assert #xs eq expected[d];rs,labels,keys:=PlacesUpToFullAutM2ThirtyFour(X,xs);K:=Parent(Evaluate(X`_coordinates[1],xs[1]));pol:=[Integers()!a:a in Coefficients(DefiningPolynomial(K,GF(3)))];ok:=[M2ThirtyFourOrbitKeys(X,x,K):x in rs];Append(~coverage,[*d,pol,keys,labels,ok*]);
 for i in [1..#rs] do ky:=MDM2FrobeniusKey([Evaluate(f,rs[i]):f in X`_coordinates],K,3,d);Append(~source,[*d,i,pol,ky*]);end for;
 if d le 5 then for j in [1..#keys] do Append(~small,[*d,j,pol,keys[j]*]);end for;end if;
 print "C234_D9_SOURCE_ORBITS",d,#xs,#rs,Cputime(tt);
end for;
Write("C2x32_C2x34_C2x40/data/C2x34/source_coverage.txt",Sprint([*source,small,coverage*]):Overwrite:=true);
g:=Genus(X);assert g eq 49;fp:=FrobeniusPolynomial(JH(68,[35]),3);assert Degree(fp) eq 2*g;tr:=[];closed:=[];
for n in [1..9] do z:=-n*Coefficient(fp,2*g-n);for j in [1..n-1] do z-:=Coefficient(fp,2*g-j)*tr[n-j];end for;Append(~tr,Integers()!z);num:=3^n+1-Integers()!z;for k in Divisors(n) do if k lt n then num-:=k*closed[k];end if;end for;assert num mod n eq 0;Append(~closed,num div n);end for;
cs:=Cusps(X);nc:=[closed[d]-#[c:c in cs|Degree(c) eq d]:d in [1..9]];assert nc eq expected;
Write("C2x32_C2x34_C2x40/data/C2x34/source_frobenius.txt",Sprint([*Coefficients(fp),[Degree(c):c in cs],nc*]):Overwrite:=true);print "C234_D9_SOURCE_ENUMERATION_DONE complete=true",Cputime(tt);quit;
