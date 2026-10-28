SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);tt:=Cputime();
load "P73/quotient19/finite_common.m";
P:=MinimalPolynomial(jAF);dj:=Degree(P);assert IsDivisibleBy(Degree(AF),dj);dd:=LCM([Denominator(c):c in Eltseq(P)]);cc:=[PF!(dd*c):c in Eltseq(P)];gg:=GCD(cc);cc:=[ExactQuotient(c,gg):c in cc];dv:=Max([Degree(c):c in cc]);assert dv gt 0 and dv le 296;assert &+[AF!c*jAF^(i-1):i->c in cc] eq 0;
out:="P73/data/quotient19/jbase_relation_data.m";
PrintFile(out,Sprintf("JBaseDegreeJ:=%o;JBaseDegreeParameter:=%o;JBasePolynomialCoefficients:=%o;",dj,dv,[[Integers()!a:a in Eltseq(c)]:c in cc]):Overwrite:=true);PrintFile(out,"print \"NONIC73_JBASE_RELATION_READY\";");printf "NONIC73_JBASE_RELATION_CERTIFIED degree_j=%o degree_parameter=%o exact_field_identity=true time=%o\n",dj,dv,Cputime(tt);quit;
