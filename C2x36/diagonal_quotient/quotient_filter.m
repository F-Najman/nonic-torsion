SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
load "C2x36/diagonal_quotient/cusp_units_model.m";
load "C2x36/diagonal_quotient/cusp_permutations.m";
Attach("common/fast_diamond_hash_m2.m");Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
function KeyThree(vs,K,p,d)
 keys:=[];for j in [1..d] do Append(~keys,&cat[[Integers()!a:a in Eltseq(v,GF(p))]:v in vs]);vs:=[v^p:v in vs];end for;return Sort(keys)[1];
end function;
function InverseThree(coords,values)
 F:=Parent(coords[1]);K:=Parent(values[1]);k:=PrimeField(K);p:=#k;e:=LCM([Degree(MinimalPolynomial(v,k)):v in values]);found:=false;
 for bits in [[1,0,0],[0,1,0],[0,0,1]] cat [[1,a,b]:a,b in [0..p-1]] cat [[0,1,a]:a in [0..p-1]] do
  gen:=&+[K!bits[i]*values[i]:i in [1..3]];pol:=MinimalPolynomial(gen,k);if Degree(pol) eq e then chosen:=bits;found:=true;break;end if;
 end for;assert found;M:=Matrix(k,[Eltseq(gen^i,k):i in [0..e-1]]);assert Rank(M) eq e;
 cc:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in values];assert &and[&+[K!cc[j][i+1]*gen^i:i in [0..e-1]] eq values[j]:j in [1..3]];
 gg:=&+[F!chosen[i]*coords[i]:i in [1..3]];fs:=[Evaluate(pol,gg)] cat [coords[j]-&+[F!cc[j][i+1]*gg^i:i in [0..e-1]]:j in [1..3]];
 O:=MaximalOrderFinite(F);I:=ideal<O|fs>;fa:=Factorization(I);places:=[Place(t[1]):t in fa|t[2] gt 0];
 if #places eq 0 then Oi:=MaximalOrderInfinite(F);Ii:=ideal<Oi|fs>;fi:=Factorization(Ii);places:=[Place(t[1]):t in fi|t[2] gt 0];end if;
 assert #places eq 1 and Degree(places[1]) eq e;assert &and[f eq 0 or Valuation(f,places[1]) gt 0:f in fs];return places[1];
end function;
function M2ThirtySixPair(E,P,Q)
 aa:=aInvariants(E);xx:=Q[1]-P[1];yy:=Q[2];return [xx^3/yy^2,(aa[2]+3*P[1])*xx^2/yy^2];
end function;
function Projection(E,P,Q)
 Q12:=3*Q;P3:=4*Q12;lambda:=Q12[2]/(Q12[1]-P[1]);bc:=M2ThirtySixPair(E,P,Q12);
 EG,iso:=IsogenyFromKernel(E,PolynomialRing(BaseRing(E)).1-P3[1]);S:=iso(4*Q);w:=(S[1]-P[1])/lambda^2;return bc cat [w];
end function;
K:=GF(5^3);places:=[p:p in Places(Z,3)|Valuation(zj,p) ge 0];assert #places eq 32;
vals:=[MDM2ToCommonField([Evaluate(f,p):f in zcoords],K):p in places];keys:=[KeyThree(v,K,5,3):v in vals];assert #Seqset(keys) eq 32;
seen:={};reps:=[];coverage:=[* *];
for i in [1..#places] do
 if keys[i] in seen then continue;end if;Append(~reps,i);b,c,w:=Explode(vals[i]);EK:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);PK:=EK![0,0];QK:=EK![b,b];ims:={};
 for bb,cc in [0,1],a in [1,5] do
  PP:=PK+6*bb*QK;QQ:=cc*PK+a*QK;lambda:=QQ[2]/(QQ[1]-PP[1]);bc:=M2ThirtySixPair(EK,PP,QQ);ww:=(w-PP[1])/lambda^2;key:=KeyThree(bc cat [ww],K,5,3);assert key in keys;Include(~ims,key);
 end for;assert IsDisjoint(seen,ims);seen join:=ims;Append(~coverage,[*i,Setseq(ims)*]);
end for;assert #seen eq 32 and #reps eq 4;print "C236_DELTA13_CUBIC_AUT_ORBITS",reps,Cputime(tt);
c0:=cs[CuspGlobalRational[1]];assert Degree(c0) eq 1;G,_,phiZ:=ClassGroup(Z:ReductionDivisor:=Divisor(c0),Proof:=true);hc:=[phiZ(Divisor(p)-Degree(p)*Divisor(c0)):p in cs];HC:=sub<G|hc>;GG,qq:=quo<G|HC>;
acs:=[hc[CuspGamma7[i]]-hc[i]+7*hc[CuspDiamond7[CuspGamma7[i]]]-7*hc[CuspDiamond7[i]]:i in [1..#cs]];assert acs[CuspGlobalRational[1]] eq G!0;AH:=sub<G|acs>;
records:=[* *];remaining:=[];
for i in [reps[1]] do
 b,c,w:=Explode(vals[i]);EK:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);PK:=EK![0,0];Q12:=EK![b,b];found:=false;
 for mult in [1,3] do
  L:=GF(5^(3*mult));EL:=EllipticCurve([L!a:a in aInvariants(EK)]);PL:=EL![0,0];Q12L:=EL![L!b,L!b];T:=4*Q12L;EP,iso:=IsogenyFromKernel(EL,PolynomialRing(L).1-T[1]);psi:=DivisionPolynomial(EL,9);
  for rt in Roots(psi) do
   x:=rt[1];aa:=aInvariants(EL);yy:=x^3+aa[2]*x^2+aa[4]*x+aa[5];if not IsSquare(yy) then continue;end if;RR:=EL![x,Sqrt(yy)];if 3*RR eq EL!0 then continue;end if;
   if (iso(RR))[1] ne L!w then continue;end if;if 3*RR eq -T then RR:=-RR;end if;assert 3*RR eq T;
   Q36:=3*(9*Q12L)+7*RR;assert 3*Q36 eq Q12L and Order(Q36) eq 36 and PL ne 18*Q36;found:=true;break;
  end for;if found then break;end if;
 end for;assert found;d:=Degree(L);assert d in [3,9];orig:=Projection(EL,PL,Q36);assert InverseThree(zcoords,orig) eq places[i];
 liftdegree:=d div 3;D:=DivisorGroup(Z)!0;images:=MDIsogeniesFastPoints(EL,[PL,Q36],7);assert &+[t[3]:t in images] eq 8;
 for t in images do
  EG,pts,e:=Explode(t);vv:=Projection(EG,pts[1],pts[2]);z:=InverseThree(zcoords,vv);n:=d*e div Degree(z);assert d*e mod Degree(z) eq 0 and n mod liftdegree eq 0;D+:=(n div liftdegree)*Divisor(z);
 end for;assert Degree(D) eq 24;
 ys:=[];ms:=[];fullmaps:=MDIsogenies(EL,7);assert #fullmaps eq 8;
 for fulliso in fullmaps do
  EG:=Codomain(fulliso);pp:=fulliso(PL);qqq:=fulliso(Q36);zz:=InverseThree(zcoords,Projection(EG,pp,qqq));jj:=Index(ys,zz);if jj eq 0 then Append(~ys,zz);Append(~ms,1);else ms[jj]+:=1;end if;
 end for;
 ref:=DivisorGroup(Z)!0;for jj in [1..#ys] do mm:=ms[jj]*3/Degree(ys[jj]);assert Denominator(mm) eq 1;ref+:=Integers()!mm*Divisor(ys[jj]);end for;assert ref eq D and Degree(ref) eq 24;print "C236_DELTA13_UNTOUCHED_FULL_NATIVE_T7_CONTROL_EQUAL degree=24",Cputime(tt);
dd:=InverseThree(zcoords,Projection(EL,PL,7*Q36));assert Degree(dd) eq 3;
 A7:=D-Divisor(places[i])-7*Divisor(dd);assert Degree(A7) eq 0;rr:=Basis(A7:Reduction:=false,Simplification:="None");principal:=#rr eq 1;assert #rr in [0,1];if principal then Append(~remaining,i);end if;
 print "C236_DELTA13_A7_CLASS_ORDERS",i,Order(qq(phiZ(Divisor(places[i])-3*Divisor(c0)))),Order(phiZ(A7));
 a:=phiZ(Divisor(places[i])-3*Divisor(c0));aa:=phiZ(A7);assert Order(qq(a)) eq 8;possible:=[n:n in [0..7]|n*aa in AH];
 print "C236_DELTA13_JACOBIAN_IMAGE_COSETS",Invariants(HC),Invariants(AH),possible;
 Ffree:=FreeAbelianGroup(#hc+1);hIn:=hom<Ffree->G|hc cat [a]>;hOut:=hom<Ffree->G|acs cat [aa]>;
 assert &*[v:v in Invariants(G)|v ne 0] eq #Image(hIn);assert &and[hOut(k) eq G!0:k in Generators(Kernel(hIn))];
 KR:=sub<G|[hIn(k):k in Generators(Kernel(hOut))]>;assert #KR*#Image(hOut) eq #Image(hIn);assert #Image(hOut) eq 8 and KR eq HC;
 Hratgens:=[&+[hc[j]:j in [1..#cs]|sigs[j] eq v]:v in Seqset(sigs)];Hrat:=sub<G|Hratgens>;print "C236_DELTA13_GLOBAL_RATIONAL_CUSP_GROUP",Invariants(Hrat),Hrat eq HC;
 print "C236_DELTA13_FULL_FINITE_A7_HOM_PASS",Invariants(KR),Invariants(Image(hOut)),"index",#Image(hOut);
 HomGroupData:=[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(h):h in hc],[Eltseq(h):h in acs],Eltseq(a),Eltseq(aa),[Eltseq(h):h in Generators(KR)],Invariants(KR),Invariants(Image(hOut)),[Eltseq(h):h in Hratgens],Invariants(Hrat)*];
 Write("C2x36/data/diagonal_quotient/quotient_filter_hom_arithmetic.txt",Sprint(HomGroupData):Overwrite:=true);

 Write("C2x36/data/diagonal_quotient/quotient_filter_first_control_arithmetic.txt",Sprint([*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(h):h in hc],[Eltseq(h):h in acs],Eltseq(a),Eltseq(aa),possible*]):Overwrite:=true);

 supp,mm:=Support(D);Append(~records,[*i,d,[*[*Degree(supp[j]),mm[j],Sprint(supp[j])*]:j in [1..#supp]*],#rr*]);
 print "C236_DELTA13_CUBIC_T7_CASE",i,d,principal,Cputime(tt);
 Write("C2x36/data/diagonal_quotient/quotient_filter_first_hecke_arithmetic.txt",Sprint([*coverage,records,remaining*]):Overwrite:=true);
end for;
/* Necessary rational-Jacobian-image filter on all2064 fresh first survivors. */
load "common/class_group_direct_RR_checked_v4_fastfb.m";
Zquery:=RationalExtensionRepresentation(Z);coordsQuery:=[Zquery!f:f in zcoords];baseQuery:=DivisorGroup(Zquery)!Divisor(c0);
Gfast,_,phiFast,witnessFilter:=MDClassGroupDirectRRCheckedV4FastFB(Zquery,"C2x36/logs/diagonal_quotient/quotient_filter_RR_identity.log");assert Gfast cmpeq G;
first:=eval Read("C2x36/data/filter_quotient18_linear_arithmetic.txt");singles:=first[3];small:=first[4];survivors:=first[7];assert #survivors eq 2064;
own:=AssociativeArray();sm:=AssociativeArray();for d in [3..9] do own[d]:=[r:r in singles|r[1] eq d];end for;for d in [3..6] do sm[d]:=[r:r in small|r[1] eq d];end for;
R5:=PolynomialRing(GF(5));fields:=AssociativeArray();for a in first[2] do fields[a[1]]:=ext<GF(5)|R5!a[2]>;end for;
cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
procedure PushProjection(~cl,~cache,~images,row)
 d:=row[1];key:=row[3];ky:=[d] cat key;if IsDefined(cache,ky) then cl:=cache[ky];return;end if;
 L:=fields[d];zz:=L.1;b:=Evaluate(R5!key[1..d],zz);c:=Evaluate(R5!key[d+1..2*d],zz);EE:=EllipticCurve([L|0,c,0,b*(1-b-c),0]);PP:=EE![0,0];QQ:=EE![b,b];assert 36*QQ eq EE!0 and 18*QQ ne EE!0 and 12*QQ ne EE!0 and PP ne 18*QQ;
 vv:=Projection(EE,PP,QQ);pp:=InverseThree(coordsQuery,vv);e:=Degree(pp);assert d mod e eq 0 and d div e in [1,3];w:=d div e;cl:=w*phiFast(Divisor(pp)-e*baseQuery);cache[ky]:=cl;
 Append(~images,[*d,key,e,w,Eltseq(cl),KeyThree(vv,L,5,d)*]);
end procedure;
for k in [1..#survivors] do
 desc:=survivors[k];rows:=[*own[desc[2]][desc[3]]*];if desc[1] ge 2 then Append(~rows,sm[desc[4]][desc[5]]);end if;if desc[1] eq 3 then Append(~rows,sm[desc[6]][desc[7]]);end if;
 cl:=G!0;for row in rows do a:=G!0;PushProjection(~a,~cache,~images,row);cl+:=a;end for;yes:=cl in HC;Append(~tests,[*k,desc,Eltseq(cl),yes*]);if yes then Append(~remaining,desc);end if;
 if k mod 50 eq 0 then
  Write("C2x36/data/diagonal_quotient/quotient_filter_arithmetic.txt",Sprint([*HomGroupData,images,tests,remaining*]):Overwrite:=true);
  print "C236_D9_GENUS9_QUOTIENT_FILTER_PROGRESS",k,#remaining,Cputime(tt);
 end if;
end for;
Write("C2x36/data/diagonal_quotient/quotient_filter_arithmetic.txt",Sprint([*HomGroupData,images,tests,remaining*]):Overwrite:=true);fprintf witnessFilter,"]\n";Flush(witnessFilter);
printf "C236_D9_GENUS9_QUOTIENT_FILTER_DONE tested=%o remaining=%o complete=true cpu=%o\n",#tests,#remaining,Cputime(tt);quit;
