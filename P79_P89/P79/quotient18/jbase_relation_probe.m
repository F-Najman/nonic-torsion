SetColumns(0);SetQuitOnError(true);SetMemoryLimit(32*10^9);tt:=Cputime();
load "P79_P89/P79/quotient18/finite_common.m";
P:=MinimalPolynomial(jAF);assert Degree(P) eq 13;dd:=LCM([Denominator(c):c in Eltseq(P)]);cc:=[PF!(dd*c):c in Eltseq(P)];gg:=GCD(cc);cc:=[ExactQuotient(c,gg):c in cc];
assert Max([Degree(c):c in cc]) eq 240;
assert &+[AF!c*jAF^(i-1):i->c in cc] eq 0;
out:="P79_P89/data/P79/quotient18/jbase_relation_data.m";
PrintFile(out,Sprintf("JBaseDegreeJ:=13;JBaseDegreeParameter:=240;JBasePolynomialCoefficients:=%o;",[[Integers()!a:a in Eltseq(c)]:c in cc]):Overwrite:=true);
PrintFile(out,"print \"NONIC79_JV_RELATION_READY\";");printf "NONIC79_JV_RELATION_CERTIFIED degree_j=13 degree_parameter=240 exact_field_identity=true time=%o\n",Cputime(tt);quit;
