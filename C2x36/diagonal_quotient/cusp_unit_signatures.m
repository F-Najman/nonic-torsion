SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();load "C2x36/diagonal_quotient/model.m";
EZ:=EllipticCurve([Z!a:a in aInvariants(E)]);PZ:=EZ![Z!a:a in Eltseq(P)];QZ:=EZ![Z!a:a in Eltseq(Q)];TZ:=4*QZ;EPZ,phi3:=IsogenyFromKernel(EZ,PolynomialRing(Z).1-TZ[1]);
aa:=aInvariants(EPZ);ok,yy:=IsSquare(root^3+aa[2]*root^2+aa[4]*root+aa[5]);assert ok;S:=EPZ![root,yy];dual:=DualIsogeny(phi3);if dual(S) eq -TZ then S:=-S;end if;assert dual(S) eq TZ;
Qmid:=phi3(9*QZ)+S;Pmid:=phi3(PZ);assert 12*Qmid eq EPZ!0 and 6*Qmid ne EPZ!0 and 4*Qmid ne EPZ!0 and Pmid ne 6*Qmid;
EPP,phi9mid:=IsogenyFromKernel(EPZ,PolynomialRing(Z).1-root);Pend:=phi9mid(Pmid);Qend:=phi9mid(phi3(9*QZ));assert 4*Qend eq EPP!0 and 2*Qend ne EPP!0 and Pend ne 2*Qend;
function Units(EE,PP,QQ)
 a:=aInvariants(EE);xx:=QQ[1]-PP[1];yy:=QQ[2];lam:=yy/xx;bb:=xx^3/yy^2;cc:=(a[2]+3*PP[1])*xx^2/yy^2;dd:=bb*(1-bb-cc);return [bb,dd,cc^2-4*dd,Discriminant(EE)/lam^12];
end function;
fs:=Units(EZ,PZ,QZ) cat [Discriminant(EPP)/Discriminant(EZ)] cat Units(EPP,Pend,Qend) cat Units(EPZ,Pmid,Qmid);
cs:=Poles(zj);sigs:=[[Valuation(f,p):f in fs]:p in cs];assert #Seqset(sigs) eq 16;
ord:=Sort([<sigs[i],Degree(cs[i]),i>:i in [1..#cs]]);cs:=[cs[t[3]]:t in ord];sigs:=[t[1]:t in ord];
Write("C2x36/data/diagonal_quotient/cusp_unit_signatures.txt",Sprint([*[Degree(p):p in cs],sigs*]):Overwrite:=true);print "C236_DELTA13_CUSP_UNIT_SIGNATURES_DONE",#cs,#Seqset(sigs),Cputime(tt);quit;
