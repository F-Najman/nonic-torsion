// Root-allocated elliptic Hecke/jet stage; no function-field normalization.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");Attach("common/fast_hecke.m");
load "C85_C91_C93/data/N93_source_units/canonical_unit_data.m";
function Code(z)return [Integers()!a:a in Eltseq(z,GF(2))];end function;
function SimplifyPoint(E,P)
 F:=BaseRing(E);a1,a2,a3,a4,a6:=Explode(aInvariants(E));assert a1 eq 1 and a3 eq 0;EE:=EllipticCurve([F!1,a2,F!0,F!0,a6+a4^2]);PP:=EE![P[1],P[2]+a4];assert jInvariant(EE)eq jInvariant(E)and 93*PP eq EE!0;assert TateNormalForm_xy(EE,PP)eq TateNormalForm_xy(E,P);return EE,PP;
end function;
function UnitValues(aa,xx,yy)
 a1,a2,a3,a4,a6:=Explode(aa);aa3:=2*yy+a3+a1*xx;aa2:=3*xx+a2;aa4:=3*xx^2+2*xx*a2+a4-a1*yy;aaa1:=2*aa4/aa3+a1;aaa2:=aa2-(aa4/aa3)^2-a1*aa4/aa3;b:=-aaa2^3/aa3^2;c:=-(aaa1*aaa2-aa3)/aa3;r:=b/c;s:=c^2/(b-c);t:=r*s-2*r+1;return[(s-r)/t,t/(s^2-s-r+1)];
end function;
function UnitMatrix(F,x0,y0)
 return VerticalJoin(Matrix(F,[[F!v[i]+x0*v[261+i]:i in[1..261]]:v in Rows(UnitXRelations)]),Matrix(F,[[F!v[i]+y0*v[261+i]:i in[1..261]]:v in Rows(UnitYRelations)]));
end function;
function CanonicalPoint(E,P)
 F:=BaseRing(E);xy:=TateNormalForm_xy(E,P);RR:=UnitMatrix(F,xy[1],xy[2]);assert Rank(RR)eq 260;zz:=Eltseq(Basis(Kernel(Transpose(RR)))[1]);ki:=Index([z ne 0:z in zz],true);return[z/zz[ki]:z in zz],xy;
end function;
function Psi93(xx,aa6)
 PS:=Parent(xx);pp:=[PS!0,PS!1,xx,xx^4+xx^3+aa6,xx^6+aa6*xx^2];
 for n in[5..93]do m:=n div 2;if IsOdd(n)then pn:=pp[m+3]*pp[m+1]^3+pp[m]*pp[m+2]^3;else pn:=pp[m+1]/xx*(pp[m+3]*pp[m]^2+pp[m-1]*pp[m+2]^2);end if;Append(~pp,pn);end for;return pp[94];
end function;
function CanonicalJets(E,P,z0)
 F:=BaseRing(E);aa:=aInvariants(E);assert aa[1]eq 1 and aa[3]eq 0 and aa[4]eq 0;PS<u>:=PowerSeriesRing(F,5);xx:=PS!P[1];yy:=PS!P[2];a2:=aa[2];a6:=PS!aa[5]+u;der:=Coefficient(Psi93(xx+u,PS!aa[5]),1);assert der ne 0 and P[1]ne 0;
 for n in[1..4]do xx+:=(Coefficient(Psi93(xx,a6),n)/der)*u^n;assert Valuation(Psi93(xx,a6))ge n+1;end for;
 for n in[1..4]do er:=Coefficient(yy^2+xx*yy+xx^3+a2*xx^2+a6,n);yy+:=(er/P[1])*u^n;assert Valuation(yy^2+xx*yy+xx^3+a2*xx^2+a6)ge n+1;end for;
 xy:=UnitValues([PS!1,PS!a2,PS!0,PS!0,a6],xx,yy);RM:=[UnitMatrix(F,Coefficient(xy[1],0),Coefficient(xy[2],0))];
 for n in[1..4]do Append(~RM,VerticalJoin(Matrix(F,[[Coefficient(xy[1],n)*v[261+i]:i in[1..261]]:v in Rows(UnitXRelations)]),Matrix(F,[[Coefficient(xy[2],n)*v[261+i]:i in[1..261]]:v in Rows(UnitYRelations)])));end for;
 ki:=Index([z ne 0:z in z0],true);assert z0[ki]eq 1;uk:=[i:i in[1..261]|i ne ki];M0:=Matrix(F,[[RM[1][r,c]:c in uk]:r in[1..Nrows(RM[1])]]);assert Rank(M0)eq 260;LL:=Solution(M0,IdentityMatrix(F,260));assert LL*M0 eq IdentityMatrix(F,260);zv:=[Vector(F,z0)];
 for n in[1..4]do
  br:=Vector(F,[0:i in[1..Nrows(M0)]]);for j in[1..n]do br+:=Vector(F,[&+[RM[j+1][r,c]*zv[n-j+1][c]:c in[1..261]]:r in[1..Nrows(M0)]]);end for;ur:=br*Transpose(LL);zz:=Vector(F,[0:i in[1..261]]);for j in[1..260]do zz[uk[j]]:=ur[j];end for;Append(~zv,zz);
 end for;
 for n in[0..4]do for r in[1..Nrows(M0)]do assert &+[&+[RM[j+1][r,c]*zv[n-j+1][c]:c in[1..261]]:j in[0..n]]eq 0;end for;end for;
 return <[Code(Coefficient(xx,n)):n in[0..4]],[Code(Coefficient(yy,n)):n in[0..4]],[[Code(zv[n+1][i]):n in[0..4]]:i in[1..261]]>;
end function;
function RecordPoint(E,P)
 E,P:=SimplifyPoint(E,P);F:=BaseRing(E);z,xy:=CanonicalPoint(E,P);jt:=CanonicalJets(E,P,z);d:=LCM([Degree(MinimalPolynomial(x,GF(2))):x in xy]);assert LCM([Degree(MinimalPolynomial(x,GF(2))):x in z])eq d;
 rec:=<[Integers()!c:c in Eltseq(DefiningPolynomial(F,GF(2)))],[Code(a):a in aInvariants(E)],[Code(P[i]):i in[1,2]],[Code(a):a in z],[Code(a):a in xy],jt>;return rec,d;
end function;
load "C85_C91_C93/data/N93_f9/records.m";
assert Complete and FieldDegree eq 9 and #Records eq 3;
PF:=PolynomialRing(GF(2));F<a>:=ext<GF(2)|PF!FieldPolynomial>;
rec:=Records[1];E:=EllipticCurve([F!PF!v:v in rec[2]]);PP:=E![F!PF!v:v in rec[3]];
assert #E eq 558 and Order(PP)eq 93 and Degree(MinimalPolynomial(jInvariant(E),GF(2)))eq 3;
SourceRecords:=[* *];SourceHeckeRecords:=[* *];SourceJPolynomials:=[* *];
sr,sd:=RecordPoint(E,PP);assert sd eq 9;Append(~SourceRecords,sr);Append(~SourceJPolynomials,[Integers()!c:c in Eltseq(MinimalPolynomial(jInvariant(E),GF(2)))]);
fast:=MDIsogeniesFast(E,PP,5);assert &+[a[3]:a in fast]eq 6;
untouched:=MDIsogenies(E,5);assert #untouched eq 6;alljs:=[jInvariant(Codomain(phi)):phi in untouched];assert #Seqset(alljs)eq 6;
trs:=[* *];
for dat in fast do
 ET,PT,e:=Explode(dat);assert Order(PT)eq 93;tr,td:=RecordPoint(ET,PT);assert (9*e)mod td eq 0;mult:=9*e div td;Append(~trs,<e,td,mult,tr>);
 printf "NONIC93_HEC5_TARGET degree=%o norm_multiplicity=%o time=%o\n",td,mult,Cputime(tt);
end for;
assert &+[r[2]*r[3]:r in trs]eq 54;Append(~SourceHeckeRecords,trs);
printf "NONIC93_ELLIPTIC_HECKE_JETS_SOURCE source_degree9=true target_closed_degrees=%o norm_multiplicities=%o all_unit_recovery_ranks260=true psi93_derivative_nonzero=true jet_precision5=true time=%o\n",[r[2]:r in trs],[r[3]:r in trs],Cputime(tt);
out:="C85_C91_C93/data/N93_source_units/elliptic_hecke5_v2_data.m";PrintFile(out,Sprintf("SourceJPolynomials:=%o;SourceRecords:=%o;SourceHeckeRecords:=%o;",SourceJPolynomials,SourceRecords,SourceHeckeRecords):Overwrite:=true);PrintFile(out,"print \"NONIC93_ELLIPTIC_HECKE_JETS_DATA_READY\";");print "NONIC93_ELLIPTIC_HECKE_JETS_STAGE2_COMPLETE";quit;
