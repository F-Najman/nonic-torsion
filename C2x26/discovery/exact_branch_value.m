/* Exact formal branch through (u,v)=(-1,0), no curve normalization. */
SetColumns(0);SetQuitOnError(true);
R<U,V>:=PolynomialRing(Rationals(),2);u:=U;v:=V;
load "mdmagma/models_X1_m_n/X1_2_26.txt";
F:=X;Ftrans:=Evaluate(F,[-1+V*U,V]);mons:=Monomials(Ftrans);coefs:=Coefficients(Ftrans);assert Minimum([Exponents(m)[2]:m in mons]) eq 3;
H:=&+[coefs[i]*U^Exponents(mons[i])[1]*V^(Exponents(mons[i])[2]-3):i in [1..#mons]];
assert Evaluate(H,[U,R!0]) eq -U^2*(U+1);dH:=Derivative(H,1);
S<z>:=LaurentSeriesRing(Rationals(),80);ww:=S!(-1);for j in [1..8] do ww-:=Evaluate(H,[ww,z])/Evaluate(dH,[ww,z]);end for;assert Valuation(Evaluate(H,[ww,z])) ge 70;
u:=-1+z*ww;v:=z;f:=eval Read("C2x26/data/discovery/miller_function_F0.txt");assert Valuation(f) eq 0;lambda:=Coefficient(f,0);print "C226_EXACT_FORMAL_PENCIL_CUSP_VALUE",lambda;
Write("C2x26/data/discovery/exact_pencil_cusp_value.txt",Sprint(lambda):Overwrite:=true);Write("C2x26/data/discovery/exact_B7_branch_series.txt",Sprint(ww):Overwrite:=true);print "C226_EXACT_FORMAL_BRANCH_VALUE_DONE complete=true";quit;
