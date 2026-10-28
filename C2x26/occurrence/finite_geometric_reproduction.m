if not assigned C226Out then C226Out:="C2x26/data/occurrence";end if;
/* Canonical self-contained existence certificate: NO finite classgroup
   construction and NO characteristic0 functionfield normalization. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
labels:=[-22,-18,-14,-10,-6,-2,1,2,3,4,5,6,7,8,9,10,11,12];
Avals:=[3,-2,4,-18,1,6,-10,11,10];Amults:=[2,2,2,1,1,1,1,1,1];Bvals:=[7,8,-22];
Cvals:=[2,5,-14,-6,9,12];Cmults:=[3,3,3,1,1,1];records:=[* *];
for p in [3,7] do
 X:=MDX11(2,26,GF(p));C:=Curve(X);assert Genus(C) eq 25;FF:=FunctionField(C);AF,am:=AlgorithmicFunctionField(FF);
 E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];b:=Q[1];assert b eq FF!X`_coordinates[1];
 rat:=[c:c in Cusps(X)|Degree(c) eq 1];assert #rat eq 18;ss:=[Valuation(b,c):c in rat];assert Sort(ss) eq labels;
 cs:=AssociativeArray();for i in [1..18] do cs[ss[i]]:=rat[i];end for;
 A:=&+[Amults[i]*Divisor(cs[Avals[i]]):i in [1..#Avals]];B:=&+[Divisor(cs[v]):v in Bvals];CC:=&+[Cmults[i]*Divisor(cs[Cvals[i]]):i in [1..#Cvals]];L:=A-B;
 assert Degree(A) eq 12 and Degree(B) eq 3 and Degree(CC) eq 12 and Degree(L) eq 9 and A ne CC;
 dimA:=Dimension(RiemannRochSpace(A));assert dimA eq 2;
 rel:=FunctionFieldDivisor(CC-A);relbasis:=Basis(rel:Reduction:=false,Simplification:="None");assert #relbasis eq 1 and Divisor(relbasis[1])+rel eq DivisorGroup(AF)!0;
 Write(C226Out cat "/cuspidal_relation_F" cat IntegerToString(p) cat ".txt",Sprint(relbasis[1]):Overwrite:=true);
 LL:=FunctionFieldDivisor(L);BB:=Basis(LL:Reduction:=false,Simplification:="None");assert #BB eq 1;
 T:=LL+Divisor(BB[1]);assert IsEffective(T) and Degree(T) eq 9;ps,ms:=Support(T);
 localdata:=[* *];
 for i in [1..#ps] do
  pp:=ps[i];ai:=[Evaluate(am(FF!f),pp):f in X`_E];EK:=EllipticCurve(ai);assert Discriminant(EK) ne 0;
  PK:=EK![Evaluate(am(FF!f),pp):f in X`_P];QK:=EK![Evaluate(am(FF!f),pp):f in X`_Q];
  assert 2*PK eq EK!0 and PK ne EK!0 and 26*QK eq EK!0 and 13*QK ne EK!0 and 2*QK ne EK!0 and PK ne 13*QK;
  K:=BaseRing(EK);assert #K eq p^Degree(pp);bc:=[Evaluate(am(FF!f),pp):f in X`_coordinates];
  pol:=[Integers()!a:a in Coefficients(DefiningPolynomial(K,GF(p)))];vals:=[[Integers()!a:a in Eltseq(z,GF(p))]:z in bc];
  Append(~localdata,[*Degree(pp),ms[i],#EK,pol,vals*]);
 end for;
 if p eq 3 then
  assert #ps eq 1 and ms eq [1] and Degree(ps[1]) eq 9 and localdata[1][3] eq 19604;
  HP:=13*Q;HQ:=P+3*Q;bh:=(HQ[1]-HP[1])^3/HQ[2]^2;
  hperm:=[Valuation(bh,cs[v]):v in labels];assert Sort(hperm) eq labels;
  for i in [1..#Avals] do assert Valuation(A,cs[hperm[Index(labels,Avals[i])]]) eq Amults[i];end for;
  assert Seqset([hperm[Index(labels,v)]:v in Bvals]) eq Seqset(Bvals);
  Append(~records,[*p,dimA,#relbasis,#BB,[Degree(x):x in ps],ms,localdata,hperm*]);
 else
  assert Sort([Degree(x):x in ps]) eq [2,2,2,3] and ms eq [1,1,1,1];
  assert #[r:r in localdata|r[1] eq 2 and r[3] eq 52] eq 3;
  Append(~records,[*p,dimA,#relbasis,#BB,[Degree(x):x in ps],ms,localdata*]);
 end if;
 Write(C226Out cat "/finite_geometric_arithmetic.txt",Sprint([*labels,Avals,Amults,Bvals,Cvals,Cmults,records*]):Overwrite:=true);
 print "C226_GEOMETRIC_FINITE_STAGE",p,"h0A",dimA,"h0L",#BB,"support",[Degree(x):x in ps],ms,Cputime(tt);
end for;
print "C226_GEOMETRIC_EXISTENCE_AND_FULL_TORSION_CERTIFICATE_DONE complete=true full_group=2x26 degree=9",Cputime(tt);quit;
