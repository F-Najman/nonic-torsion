SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
root:="C78_C117/data/C78_known39/";
load "C78_C117/data/C78_known39/inputs/known39.m";
R<T>:=PolynomialRing(Rationals());f:=R!FieldCoefficients;xp:=R!XCoefficients;yp:=R!YCoefficients;
assert Degree(f)eq 9 and IsIrreducible(f);K<a>:=NumberField(f);xx:=Evaluate(xp,a);yy:=Evaluate(yp,a);
function FromXY(x,y)
 r:=(x^2*y-x*y+y-1)/(x^2*y-x);s:=(x*y-y+1)/(x*y);b:=r*s*(r-1);c:=s*(r-1);return EllipticCurve([1-c,-b,-b,0,0]);
end function;
function TateBC(E,P)
 x0:=P[1];y0:=P[2];a1,a2,a3,a4,a6:=Explode(aInvariants(E));aa3:=2*y0+a3+a1*x0;aa2:=3*x0+a2;aa4:=3*x0^2+2*x0*a2+a4-a1*y0;
 aaa1:=2*aa4/aa3+a1;aaa2:=aa2-(aa4/aa3)^2-a1*aa4/aa3;b:=-aaa2^3/aa3^2;c:=-(aaa1*aaa2-aa3)/aa3;return [b,c];
end function;
function NFCode(v)return [[Numerator(z),Denominator(z)]:z in Eltseq(K!v)];end function;
function QCode(v)return [[Numerator(z),Denominator(z)]:z in v];end function;
function CodeElt(v)return [Integers()!c:c in Eltseq(v,PrimeField(Parent(v)))];end function;
function Less(a,b)
 for i in[1..#a]do if a[i]lt b[i]then return true;elif a[i]gt b[i]then return false;end if;end for;return false;
end function;
function FKey(v)
 ss:=[&cat[CodeElt(x^(2^i)):x in v]:i in[0..8]];ans:=ss[1];for w in ss do if Less(w,ans)then ans:=w;end if;end for;return ans;
end function;
function EvalMod(g,z)
 L:=Parent(z);return &+[L!Numerator(Coefficient(g,i))/(L!Denominator(Coefficient(g,i)))*z^i:i in[0..Degree(g)]];
end function;
E:=FromXY(xx,yy);P:=E![0,0];assert Discriminant(E)ne 0 and 39*P eq E!0 and 13*P ne E!0 and 3*P ne E!0;assert Degree(MinimalPolynomial(jInvariant(E)))eq 3;
assert Discriminant(f)eq -3^14*31^3;assert Numerator(Discriminant(f))mod 2 ne 0 and Numerator(Discriminant(f))mod 29 ne 0;
f2:=Polynomial(GF(2),f);assert IsIrreducible(f2);L2<z2>:=ext<GF(2)|f2>;E2:=FromXY(EvalMod(xp,z2),EvalMod(yp,z2));P2:=E2![0,0];assert Discriminant(E2)ne 0 and 39*P2 eq E2!0 and 13*P2 ne E2!0 and 3*P2 ne E2!0;
E29:=FromXY(EvalMod(xp,GF(29)!10),EvalMod(yp,GF(29)!10));P29:=E29![0,0];assert Discriminant(E29)eq 1 and #E29 eq 39 and 39*P29 eq E29!0 and 13*P29 ne E29!0 and 3*P29 ne E29!0;
units:=[i:i in[1..19]|GCD(i,39)eq 1];assert #units eq 12;GB:=[TateBC(E,m*P):m in units];FB:=[TateBC(E2,m*P2):m in units];keys:=[FKey(v):v in FB];
for i in[1..12]do
 assert GB[i]eq TateBC(E,(39-units[i])*P);
 assert &and[Denominator(c)mod 2 ne 0:c in &cat[Eltseq(v):v in GB[i]]];red:=[EvalMod(R!Eltseq(v),z2):v in GB[i]];assert red eq FB[i];
end for;
AA:=[r[1]:r in Roots(Polynomial(K,f))];autos:=[hom<K->K|r>:r in AA];assert #autos eq 3;
AE:=[];for rt in AA do assert &and[Denominator(c)mod 2 ne 0:c in Eltseq(rt)];rr:=EvalMod(R!Eltseq(rt),z2);inds:=[i:i in[0..8]|rr eq z2^(2^i)];assert #inds eq 1;Append(~AE,inds[1]);end for;assert Sort(AE)eq[0,3,6];
globalgroups:=[];
for i in[1..12]do
 group:=[m:m in[1..12]|exists{phi:phi in autos|[phi(v):v in GB[i]]eq GB[m]}];finite:=[m:m in[1..12]|keys[m]eq keys[i]];assert group eq finite and #group eq 3;Append(~globalgroups,[units[m]:m in group]);
end for;assert #Seqset(globalgroups)eq 4 and #Seqset(keys)eq 4;
PrintFile(root cat"global39.json",Sprintf("{\"complete\":true,\"point_order\":39,\"field_degree\":9,\"j_degree\":3,\"polynomial_discriminant\":%o,\"global_A_invariants\":%o,\"global_point\":[0,0],\"field_polynomial\":%o,\"full_torsion_prime\":29,\"residue_root\":10,\"finite29_a_invariants\":%o,\"finite29_discriminant\":1,\"finite29_order\":39,\"finite2_order\":%o,\"global_automorphism_generators\":%o,\"automorphism_Frobenius_exponents\":%o,\"physical_multipliers\":%o,\"global_trace_groups\":%o,\"finite_closed_orbits\":4}",Discriminant(f),[NFCode(v):v in aInvariants(E)],QCode(Coefficients(f)),[Integers()!c:c in aInvariants(E29)],#E2,[NFCode(v):v in AA],AE,units,Setseq(Seqset(globalgroups))):Overwrite:=true);
out:=Open(root cat"diamond39.jsonseq","w");
for i in[1..12]do
 found:=false;lambda:=0;for lam in[0..20]do pp:=MinimalPolynomial(GB[i][1]+lam*GB[i][2]);if Degree(pp)eq 9 then found:=true;lambda:=lam;break;end if;end for;assert found;
 fprintf out,"{\"multiplier\":%o,\"global_Tate_bc\":%o,\"primitive_lambda\":%o,\"primitive_polynomial\":%o,\"finite2_Tate_bc\":%o,\"finite2_frobenius_key\":%o,\"trace_group\":%o}\n",units[i],[NFCode(v):v in GB[i]],lambda,QCode(Coefficients(pp)),[CodeElt(v):v in FB[i]],keys[i],globalgroups[i];Flush(out);
end for;delete out;
printf "NONIC78_KNOWN39_GLOBAL_DIAMOND_COMPLETE degree9_exact39=true fulltors39_F29=true physical12_global4_finite4=true all_global_reductions_bound=true time=%o\n",Cputime(tt);quit;
