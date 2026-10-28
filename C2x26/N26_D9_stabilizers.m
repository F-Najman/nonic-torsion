SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_diamond_hash_m2.m");
load "C2x26/data/N26_D9_threequotients_survivors.m";
R:=PolynomialRing(GF(3));records:=[* *];
for rec in N26D9ThreeQuotientSurvivors do
 idx,pol,vals:=Explode(rec);K<z>:=ext<GF(3)|R!pol>;bc:=[Evaluate(R!v,z):v in vals];b,c:=Explode(bc);
 E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P:=E![0,0];Q:=E![b,b];twos:=[P,13*Q,P+13*Q];
 original:=MDM2FrobeniusKey(bc,K,3,9);stabs:=[* *];keys:=[];
 for pi in [1..3] do for ri in [1..3] do if pi eq ri then continue;end if;
  for a in [1..6] do
   U:=twos[pi];V:=twos[ri]+14*a*Q;r:=U[1];xx:=V[1]-r;yy:=V[2];
   vv:=[xx^3/yy^2,(c+3*r)*xx^2/yy^2];key:=MDM2FrobeniusKey(vv,K,3,9);Append(~keys,key);
   if key eq original then
    phases:=[i:i in [0..8]|vv eq [v^(3^i):v in bc]];assert #phases eq 1;
    pP:=pi in [1,3] select 1 else 0;pQ:=pi in [2,3] select 13 else 0;
    qP:=ri in [1,3] select 1 else 0;qQ:=((ri in [2,3] select 13 else 0)+14*a) mod 26;
    assert U eq pP*P+pQ*Q and V eq qP*P+qQ*Q;
    U2:=pP*U+pQ*V;V2:=qP*U+qQ*V;U3:=pP*U2+pQ*V2;V3:=qP*U2+qQ*V2;
    assert U3 eq P and V3 in [Q,-Q];
    Append(~stabs,[*pi,ri,a,pP,pQ,qP,qQ,phases[1],V3 eq Q select 1 else -1*]);
   end if;
  end for;
 end for;end for;
 assert #stabs eq 3 and #Seqset(keys) eq 12;
 Append(~records,[*idx,pol,vals,stabs,keys*]);print "N26_D9_EXACT_STABILIZER",idx,stabs,Cputime(tt);
end for;
Write("C2x26/data/N26_D9_stabilizers_arithmetic.txt",Sprint(records):Overwrite:=true);
print "N26_D9_STABILIZERS_DONE complete=true",Cputime(tt);quit;
