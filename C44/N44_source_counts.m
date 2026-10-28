SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();S:=CuspidalSubspace(ModularSymbolsH(44,[1],2,1));assert Dimension(S) eq 36;
T3:=HeckeOperator(S,3);D3:=(T3^2-HeckeOperator(S,9))/3;assert D3^10 eq IdentityMatrix(Rationals(),36);
assert Integers()!Determinant(IdentityMatrix(Rationals(),36)-T3+3*D3) eq 10863413200701600000;
M0:=2*IdentityMatrix(Rationals(),36);M1:=T3;closed:=[];
for r in [1..9] do if r eq 1 then Mr:=M1;else Mr:=T3*M1-3*D3*M0;M0:=M1;M1:=Mr;end if;nr:=Integers()!(3^r+1-Trace(Mr));cr:=(nr-&+[Integers()|dd*closed[dd]:dd in Divisors(r)|dd lt r]) div r;Append(~closed,cr);end for;
nc:=closed;nc[1]-:=15;nc[2]-:=5;nc[5]-:=3;assert nc eq [0,0,0,0,60,160,240,920,2260];
PrintFile("C44/data/N44_frobenius3_source.json",Sprintf("{\"level\":44,\"p\":3,\"genus\":36,\"local_J_order\":10863413200701600000,\"closed_place_counts\":%o,\"noncuspidal_closed_counts\":%o,\"T3_entries\":%o,\"Diamond3_entries\":%o}",closed,nc,[[Numerator(a),Denominator(a)]:a in Eltseq(T3)],[[Numerator(a),Denominator(a)]:a in Eltseq(D3)]):Overwrite:=true);
printf "NONIC44_SOURCE_FROBENIUS_COUNTS_PASS noncuspidal=%o order_correct=true time=%o\n",nc,Cputime(tt);quit;
