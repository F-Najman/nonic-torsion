assert assigned C226Out;
/* Standalone explicit degree9 torsion occurrence, no modular-curve dependency. */
SetColumns(0);SetQuitOnError(true);tt:=Cputime();
a:=eval Read("C2x26/data/explicit/example.txt");R<t>:=PolynomialRing(Rationals());f:=R!a[1];assert Degree(f) eq 9 and IsIrreducible(PolynomialRing(GF(3))!f);K<z>:=NumberField(f);
b:=Evaluate(R!a[3][1],z);c:=Evaluate(R!a[3][2],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Discriminant(E) ne 0 and 2*P eq E!0 and P ne E!0 and 26*Q eq E!0 and 13*Q ne E!0 and 2*Q ne E!0 and P ne 13*Q;
R3<t3>:=PolynomialRing(GF(3));L3<w3>:=ext<GF(3)|R3!f>;b3:=Evaluate(R3!a[3][1],w3);c3:=Evaluate(R3!a[3][2],w3);E3:=EllipticCurve([L3|0,c3,0,b3*(1-b3-c3),0]);assert Discriminant(E3) ne 0 and #E3 eq 19604;
R7<t7>:=PolynomialRing(GF(7));f7:=R7!f;assert GCD(f7,Derivative(f7)) eq R7!1;fac:=Factorization(f7);assert Sort([Degree(h[1]):h in fac]) eq [2,2,2,3];bounds:=[];
for h in fac do L<w>:=ext<GF(7)|h[1]>;bb:=Evaluate(R7!a[3][1],w);cc:=Evaluate(R7!a[3][2],w);EE:=EllipticCurve([L|0,cc,0,bb*(1-bb-cc),0]);assert Discriminant(EE) ne 0;Append(~bounds,[Degree(L),#EE]);assert #EE eq (Degree(L) eq 2 select 52 else 364);end for;
roots:=[K!0,(13*Q)[1],(P+13*Q)[1]];checks:=[];
for i in [1..3] do r:=roots[i];other:=[roots[j]:j in [1..3]|j ne i];A:=c+3*r;B:=(r-other[1])*(r-other[2]);delta:=other[1]-other[2];assert delta^2 eq A^2-4*B;values:=[2*x*(delta^2-A*x):x in [delta,-delta]];squares:=[IsSquare(v):v in values];assert squares eq [false,false];Append(~checks,[*i,[Eltseq(v):v in values],squares*]);end for;
Write(C226Out cat "/native_arithmetic.txt",Sprint([*a,bounds,checks*]):Overwrite:=true);print "C226_STANDALONE_EXPLICIT_OCCURRENCE_PASS degree=9 full_group=2x26 no4on3isogenytargets=true complete=true",Cputime(tt);quit;
