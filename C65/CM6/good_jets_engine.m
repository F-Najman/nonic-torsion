SetColumns(0);SetQuitOnError(true);SetMemoryLimit((CaseMode eq 6 select 12 else 4)*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "C65/data/source_units/canonical_unit_data.m";

load "C65/CM6/explicit3_isogenies.m";
g:=121;JPrec:=16;
function Code(z)return [Integers()!v:v in Eltseq(z,GF(2))];end function;
function SimplifyPoint(E,P)
 aa:=aInvariants(E);F:=BaseRing(E);
 if aa[1] eq 0 then assert aa[3] ne 0;return E,P;end if;
 assert aa[1] eq 1 and aa[3] eq 0;EE:=EllipticCurve([F!1,aa[2],F!0,F!0,aa[5]+aa[4]^2]);PP:=EE![P[1],P[2]+aa[4]];assert Order(PP) eq 65 and TateNormalForm_xy(EE,PP) eq TateNormalForm_xy(E,P);return EE,PP;
end function;
function UnitValues(aa,x,y)
 a1,a2,a3,a4,a6:=Explode(aa);aa3:=2*y+a3+a1*x;aa2:=3*x+a2;aa4:=3*x*x+2*x*a2+a4-a1*y;aaa1:=2*aa4/aa3+a1;aaa2:=aa2-(aa4/aa3)^2-a1*aa4/aa3;b:=-aaa2^3/aa3^2;c:=-(aaa1*aaa2-aa3)/aa3;r:=b/c;s:=c*c/(b-c);t:=r*s-2*r+1;return[(s-r)/t,t/(s*s-s-r+1)];
end function;
function UnitMatrix(F,x,y)
 return VerticalJoin(Matrix(F,[[F!v[i]+x*v[g+i]:i in[1..g]]:v in Rows(UnitXRelations)]),Matrix(F,[[F!v[i]+y*v[g+i]:i in[1..g]]:v in Rows(UnitYRelations)]));
end function;
function Psi65(x,aa)
 a1,a2,a3,a4,a6:=Explode(aa);b2:=a1*a1+4*a2;b4:=a1*a3+2*a4;b6:=a3*a3+4*a6;b8:=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4;s2:=a1*x+a3;
 s3:=3*x^4+b2*x^3+3*b4*x*x+3*b6*x+b8;s4:=s2*(2*x^6+b2*x^5+5*b4*x^4+10*b6*x^3+10*b8*x*x+(b2*b8-b4*b6)*x+b4*b8-b6*b6);pp:=[Parent(x)!0,Parent(x)!1,s2,s3,s4];
 for n in[5..65]do m:=n div 2;if IsOdd(n)then pn:=pp[m+3]*pp[m+1]^3-pp[m]*pp[m+2]^3;else pn:=pp[m+1]/s2*(pp[m+3]*pp[m]^2-pp[m-1]*pp[m+2]^2);end if;Append(~pp,pn);end for;return pp[66];
end function;
function RecordPoint(E,P)
 E,P:=SimplifyPoint(E,P);F:=BaseRing(E);aa:=aInvariants(E);xy:=TateNormalForm_xy(E,P);MM:=UnitMatrix(F,xy[1],xy[2]);assert Rank(MM) eq g-1;z0:=Eltseq(Basis(Kernel(Transpose(MM)))[1]);ki:=Index([z ne 0:z in z0],true);z0:=[z/z0[ki]:z in z0];
 PS<u>:=PowerSeriesRing(F,JPrec);aas:=[PS!v:v in aa];mode:=aa[1] eq 0 select 1 else 6;aas[mode eq 6 select 5 else 1]+:=u;xx:=PS!P[1];yy:=PS!P[2];der:=Coefficient(Psi65(xx+u,[PS!v:v in aa]),1);assert der ne 0 and Coefficient(Psi65(xx,[PS!v:v in aa]),0) eq 0;
 for n in[1..JPrec-1]do xx+:=(Coefficient(Psi65(xx,aas),n)/der)*u^n;assert Valuation(Psi65(xx,aas)) ge n+1;end for;
 dy:=aa[1]*P[1]+aa[3];assert dy ne 0;
 for n in[1..JPrec-1]do er:=yy^2+aas[1]*xx*yy+aas[3]*yy-xx^3-aas[2]*xx^2-aas[4]*xx-aas[5];yy+:=(Coefficient(er,n)/dy)*u^n;assert Valuation(yy^2+aas[1]*xx*yy+aas[3]*yy-xx^3-aas[2]*xx^2-aas[4]*xx-aas[5]) ge n+1;end for;
 uxy:=UnitValues(aas,xx,yy);RM:=[MM];
 for n in[1..JPrec-1]do Append(~RM,VerticalJoin(Matrix(F,[[Coefficient(uxy[1],n)*v[g+i]:i in[1..g]]:v in Rows(UnitXRelations)]),Matrix(F,[[Coefficient(uxy[2],n)*v[g+i]:i in[1..g]]:v in Rows(UnitYRelations)])));end for;
 uk:=[i:i in[1..g]|i ne ki];M0:=Matrix(F,[[MM[r,c]:c in uk]:r in[1..Nrows(MM)]]);LL:=Solution(M0,IdentityMatrix(F,g-1));assert LL*M0 eq IdentityMatrix(F,g-1);zv:=[Vector(F,z0)];
 for n in[1..JPrec-1]do
  br:=Vector(F,[0:i in[1..Nrows(M0)]]);for j in[1..n]do br+:=Vector(F,[&+[RM[j+1][r,c]*zv[n-j+1][c]:c in[1..g]]:r in[1..Nrows(M0)]]);end for;ur:=br*Transpose(LL);zz:=Vector(F,[0:i in[1..g]]);for j in[1..g-1]do zz[uk[j]]:=ur[j];end for;Append(~zv,zz);
 end for;
 for n in[0..JPrec-1],r in[1..Nrows(M0)]do assert &+[&+[RM[j+1][r,c]*zv[n-j+1][c]:c in[1..g]]:j in[0..n]] eq 0;end for;
 assert exists{v:v in Eltseq(zv[2])|v ne 0};degree:=LCM([Degree(MinimalPolynomial(v,GF(2))):v in xy]);assert LCM([Degree(MinimalPolynomial(v,GF(2))):v in z0]) eq degree;
 rec:=[*[Integers()!v:v in Eltseq(DefiningPolynomial(F,GF(2)))],[Code(v):v in aa],[Code(P[1]),Code(P[2])],[Code(v):v in z0],[Code(v):v in xy],mode,[[Code(Coefficient(v,n)):n in[0..JPrec-1]]:v in[xx,yy]],[[Code(zv[n+1][i]):n in[0..JPrec-1]]:i in[1..g]]*];
 return rec,degree;
end function;
if CaseMode eq 6 then inds:=[1..#Records];else assert CaseMode eq 8;inds:=[1];end if;
F:=ext<GF(2)|PolynomialRing(GF(2))!FieldPolynomial>;ev:=func<v|Evaluate(PolynomialRing(GF(2))!v,F.1)>;Sources:=[* *];Targets:=[* *];Kernels:=[* *];
for ii in inds do
 rr:=Records[ii];E:=EllipticCurve([ev(v):v in rr[2]]);P:=E![ev(v):v in rr[3]];assert Order(P) eq 65;E,P:=SimplifyPoint(E,P);ss:=[* *];
 for k in[1,3,9]do rec,dd:=RecordPoint(E,k*P);assert dd eq CaseMode;Append(~ss,rec);end for;Append(~Sources,ss);
 fast:=C65Isogs(E,P,3);trs:=[* *];krs:=[* *];assert &+[v[3]:v in fast] eq 4;
 for dat in fast do
  EG,PG,e,kpol,EL,PL:=Explode(dat);tr,td:=RecordPoint(EG,PG);tr5,td5:=RecordPoint(EG,3*PG);assert td eq td5 and (CaseMode*e) mod td eq 0;mult:=CaseMode*e div td;Append(~trs,[*e,td,mult,tr,tr5*]);
  LF:=BaseRing(EL);Append(~krs,[*[Integers()!v:v in Eltseq(DefiningPolynomial(LF,GF(2)))],[Code(v):v in aInvariants(EL)],[Code(PL[1]),Code(PL[2])],[Code(v):v in Eltseq(kpol)],[Code(v):v in aInvariants(EG)],[Code(PG[1]),Code(PG[2])],e*]);
 end for;assert &+[r[2]*r[3]:r in trs] eq CaseMode*4;Append(~Targets,trs);Append(~Kernels,krs);
 printf "NONIC65_CM6_GOOD_JETS_SOURCE case=%o index=%o target_degrees=%o norm_weights=%o unit_rank251=true genuine_tangent=true precision16=true time=%o\n",CaseMode,ii,[v[2]:v in trs],[v[3]:v in trs],Cputime(tt);
end for;
out:=Sprintf("C65/data/CM6/good%o_jets_data.txt",CaseMode);Write(out,Sprint([*CaseMode,JPrec,Sources,Targets,Kernels*]):Overwrite:=true);printf "NONIC65_CM6_GOOD_JETS_COMPLETE case=%o time=%o\n",CaseMode,Cputime(tt);quit;
