SetColumns(0);SetQuitOnError(true);SetMemoryLimit(32*10^9);tt:=Cputime();
load "P79_P89/P89/quotient27/finite_common.m";
P:=MinimalPolynomial(jAF);assert Degree(P) eq 25;dd:=LCM([Denominator(c):c in Eltseq(P)]);cc:=[PF!(dd*c):c in Eltseq(P)];gg:=GCD(cc);cc:=[ExactQuotient(c,gg):c in cc];
assert Max([Degree(c):c in cc]) eq 360;
assert &+[AF!c*jAF^(i-1):i->c in cc] eq 0;
out:="P79_P89/data/P89/quotient27/jv_relation_data.m";
PrintFile(out,Sprintf("JVDegreeJ:=25;JVDegreeV:=360;JVPolynomialCoefficients:=%o;",[[Integers()!a:a in Eltseq(c)]:c in cc]):Overwrite:=true);
PrintFile(out,"print \"NONIC89_JV_RELATION_READY\";");printf "NONIC89_JV_RELATION_CERTIFIED degree_j=25 degree_v=360 exact_field_identity=true time=%o\n",Cputime(tt);quit;
