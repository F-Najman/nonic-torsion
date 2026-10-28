/* Rational elimination from the explicit Miller pencil; no normalization. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
R<U,V>:=PolynomialRing(Rationals(),2);u:=U;v:=V;
load "mdmagma/models_X1_m_n/X1_2_26.txt";
plane:=X;F2:=FieldOfFractions(R);u:=F2!U;v:=F2!V;f:=eval Read("C2x26/data/discovery/miller_function_F0.txt");lam:=eval Read("C2x26/data/discovery/exact_pencil_cusp_value.txt");
num:=Numerator(f-lam);den:=Denominator(f-lam);print "C226_ELIMINATION_POLYNOMIALS_READY",TotalDegree(num),Degree(num,1),Degree(num,2),Cputime(tt);
res:=Resultant(plane,num,1);assert res ne 0;print "C226_ELIMINATION_RESULTANT_READY",Degree(res,2),Cputime(tt);
T<t>:=PolynomialRing(Rationals());rr:=Evaluate(res,[T!0,t]);fac:=Factorization(rr);print "C226_RESULTANT_FACTOR_DEGREES",[[Degree(a[1]),a[2]]:a in fac],Cputime(tt);
Write("C2x26/data/discovery/resultant_factors.txt",Sprint(fac):Overwrite:=true);ffs:=[a[1]:a in fac|Degree(a[1]) eq 9];assert #ffs eq 1;pol:=ffs[1];K<z>:=NumberField(pol);SK<s>:=PolynomialRing(K);
pp:=Evaluate(plane,[s,SK!z]);nn:=Evaluate(num,[s,SK!z]);g:=GCD(pp,nn);assert Degree(g) eq 1;uk:=-Coefficient(g,0)/Coefficient(g,1);assert Evaluate(den,[uk,z]) ne 0;
u:=uk;v:=z;
q:=(-u-v-2)/(u-v);t:=(u^3+u^2*v+2*u^2-u*v^2-v^3-2*v^2)/(u^3+u^2*v+2*u^2+u*v^2+2*u*v+v^3+2*v^2-2);
b:=(t+1)*(q*t+1)/t^2;c:=(t^2-2*q*t-2)/t^2;E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P2:=E![0,0];Q26:=E![b,b];assert 2*P2 eq E!0 and P2 ne E!0 and 26*Q26 eq E!0 and 13*Q26 ne E!0 and 2*Q26 ne E!0 and P2 ne 13*Q26;
print "C226_EXPLICIT_K9_AND_2x26_POINT_PASS",Cputime(tt);
Write("C2x26/data/discovery/explicit_K9_example.txt",Sprint([*Coefficients(pol),Eltseq(uk),[Eltseq(b),Eltseq(c)]*]):Overwrite:=true);
roots:=[K!0,(13*Q26)[1],(P2+13*Q26)[1]];checks:=[];
for i in [1..3] do
 r:=roots[i];others:=[roots[j]:j in [1..3]|j ne i];A:=c+3*r;B:=(r-others[1])*(r-others[2]);delta:=others[1]-others[2];assert delta^2 eq A^2-4*B;
 vals:=[2*x*(delta^2-A*x):x in [delta,-delta]];sq:=[IsSquare(v):v in vals];Append(~checks,[IsSquare(B)] cat sq);print "C226_TWO_ISOGENY_FOUR_DIVISIBILITY",i,checks[#checks];
 if &or sq then
  j:=Index(sq,true);x:=[delta,-delta][j];EP:=EllipticCurve([K|0,-2*A,0,delta^2,0]);R4:=EP![x,Sqrt(vals[j])];assert 4*R4 eq EP!0 and 2*R4 ne EP!0;
  target,ph:=IsogenyFromKernel(E,PolynomialRing(K).1-r);U:=ph(Q26);assert 26*U eq target!0 and 2*U ne target!0;
  print "C226_POSSIBLE_C52_TARGET_FOUND",i;
 end if;
end for;
Write("C2x26/data/discovery/two_isogeny_four_checks.txt",Sprint(checks):Overwrite:=true);print "C226_EXPLICIT_FIELD_AND_ISOGENIES_DONE complete=true",Cputime(tt);quit;
