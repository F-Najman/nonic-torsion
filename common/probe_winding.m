// Numeric1 denotes nonzero winding, not an L-ratio.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
if Type(N)eq MonStgElt then N:=StringToInteger(N);end if;assert N in[57,65,85,91,93,117];
S:=CuspidalSubspace(ModularSymbolsH(N,[1],2,1));us:=[u:u in[1..N]|GCD(u,N)eq 1];HH:=Seqset(us);rr:=[* *];z:=0;
printf "NONIC_CHARACTER_PROBE_START N=%o genus=%o time=%o\n",N,Dimension(S),Cputime(tt);
for A in NewformDecomposition(S)do
 chi:=DirichletCharacter(A);W:=WindingElementProjection(A);L:=IsZero(W)select 0 else 1;ker:=[u:u in us|chi(u)eq 1];assert -1 mod N in ker;
 if L eq 0 then HH meet:=Seqset(ker);z+:=1;end if;
 Append(~rr,<Dimension(A),Level(AssociatedNewSpace(A)),L,Order(chi),ker>);
 printf "NONIC_CHARACTER_FACTOR N=%o dim=%o newlevel=%o winding_nonzero=%o character_order=%o kernel_size=%o time=%o\n",N,Dimension(A),Level(AssociatedNewSpace(A)),L,Order(chi),#ker,Cputime(tt);
end for;
out:=(N in[85,91,93] select "C85_C91_C93" else N eq 117 select "C78_C117" else "C"cat IntegerToString(N))cat"/data/character_probes/N"cat IntegerToString(N)cat"_winding_data.m";
PrintFile(out,Sprintf("N:=%o;SourceGenus:=%o;ConstituentWindingFlags:=%o;ZeroLFactors:=%o;FiniteDiamondUnits:=%o;FiniteProjectiveDiamondOrder:=%o;Complete:=true;",N,Dimension(S),rr,z,Sort(Setseq(HH)),#HH div 2):Overwrite:=true);
printf "NONIC_CHARACTER_PROBE_COMPLETE N=%o genus=%o zero_L_factors=%o finite_diamond_order=%o units=%o time=%o\n",N,Dimension(S),z,#HH div 2,Sort(Setseq(HH)),Cputime(tt);quit;
