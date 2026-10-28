from pathlib import Path
import ast,json,re,sys,time
from collections import Counter
from sage.all import GF,QQ,ZZ,PolynomialRing,EllipticCurve,divisors,matrix,vector
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
S=Path(__file__).resolve().parents[1]/'data';P=S/'diagonal_genus17_lattice';t0=time.monotonic()
def rd(p):
 s=p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False')
 s=re.sub(r'\([^()\n]*F\.1[^()\n]*\)',lambda m:json.dumps(m.group()),s)
 return ast.literal_eval(s)
class Quotient:
 def __init__(self,R,H):
  self.n=len(R[0])-1;self.M=hermite_normal_form(Matrix([r[:-1] for r in R+H]).T);self.index=abs(int(self.M.det()))
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0;w=list(v[:-1])
  for j in range(self.n-1,-1,-1):
   q=w[j]//int(self.M[j,j])
   for i in range(j+1):w[i]-=q*int(self.M[i,j])
  return tuple(w)
def add(*vs):return [sum(x)for x in zip(*vs)]
def mul(n,v):return [n*x for x in v]
def projection(E,P2,Q):
 S7=6*Q;T3=14*Q;R2=21*Q;Z2=P2+R2;x,y=S7[0],S7[1];m=(3*x*x+2*E.a2()*x+E.a4())/(2*y);a2p=3*x+E.a2()-m*m;a3p=2*y
 k=a3p/a2p;b=-a2p**3/a3p**2;c=1-2*m*a2p/a3p;v=b/c;assert b==v*v*(v-1)and c==v*(v-1);H=v*v-v+1;assert H
 u=4*(P2[0]-R2[0])*(P2[0]-Z2[0])*(R2[0]-Z2[0])/(k**6*v**4*(v-1)**4);x3=(T3[0]+E.a2()/3)/(k*k*H*H);y3=T3[1]/(k**3*H**3)
 t=u*u;D=t*t+13*t+49;N=t*t+5*t+1;M=t**4+14*t**3+63*t*t+70*t-7
 assert y3*y3==x3**3-N/(48*D)*x3+M/(864*D*D)
 return [u,x3,y3]
def degree(vs,d):return min(int(e)for e in divisors(d)if all(v**(5**int(e))==v for v in vs))
R5=PolynomialRing(GF(5),'z');fcache={}
def fld(pol):
 k=tuple(pol)
 if k not in fcache:fcache[k]=GF(5**(len(pol)-1),'z',modulus=R5(pol))
 return fcache[k]
def ev(K,cs):return sum(K(c)*K.gen()**i for i,c in enumerate(cs))
def curve(d):
 pol,ai,pv,qv=d;K=fld(pol);E=EllipticCurve(K,[ev(K,a)for a in ai]);pp=[ev(K,a)for a in pv];qq=[ev(K,a)for a in qv];assert pp[2]==qq[2]==1;return K,E,E(pp),E(qq)

from sage.all import FunctionField
zdata,udata=rd(P/'parameter_data.txt')
def evalrat(cs,u):return sum(u.parent()(a)*u**i for i,a in enumerate(cs[0]))/sum(u.parent()(a)*u**i for i,a in enumerate(cs[1]))
def project17(E,P2,Q):
 vals=projection(E,P2,Q);S7=6*Q;x,y=S7[0],S7[1];m=(3*x*x+2*E.a2()*x+E.a4())/(2*y);a2p=3*x+E.a2()-m*m;a3p=2*y;kv=a3p/a2p;b=-a2p**3/a3p**2;c=1-2*m*a2p/a3p;v=b/c;H=v*v-v+1;vs=[v,(v-1)/v,1/(1-v)];rs=[(T[0]+E.a2()/3)/(kv*kv*H*H)for T in[P2,21*Q,P2+21*Q]];w=sum(r*v for r,v in zip(rs,vs));u=vals[0];t=u*u;D=t*t+13*t+49;N=t*t+5*t+1;M=t**4+14*t**3+63*t*t+70*t-7
 assert w**3-N/48*w+M*(2*t+13)/(1728*D)-u/(8*D)==0
 zz=sum(evalrat(cs,u)*w**i for i,cs in enumerate(zdata));assert evalrat(udata,zz)==u
 return[zz,vals[1],vals[2]]

from sage.all import lcm,PowerSeriesRing,identity_matrix
RU=PolynomialRing(GF(5),'u');KU=RU.fraction_field();U0=KU.gen();T0=U0*U0;D0=T0*T0+13*T0+49;N0=T0*T0+5*T0+1;M0=T0**4+14*T0**3+63*T0*T0+70*T0-7;RW=PolynomialRing(KU,'w');W0=RW.gen();fw=W0**3-N0/48*W0+M0*(2*T0+13)/(1728*D0)-U0/(8*D0);CQ=RW.quotient(fw,'w');WW=CQ.gen();Z0=sum(CQ(evalrat(c,U0))*WW**i for i,c in enumerate(zdata));un=sum(CQ(c)*Z0**i for i,c in enumerate(udata[0]));ud=sum(CQ(c)*Z0**i for i,c in enumerate(udata[1]));assert un==CQ(U0)*ud and ud!=0
assert len(udata[0])-1==3 and len(udata[1])-1==2
print('C242_GENUS17_INDEPENDENT_RATIONAL_PARAMETER_IDENTITY_PASS',flush=True)
Rz=PolynomialRing(GF(5),'z');Kz=Rz.fraction_field();z=Kz.gen();u=evalrat(udata,z);t=u*u;DD=t*t+13*t+49;NN=t*t+5*t+1;MM=t**4+14*t**3+63*t*t+70*t-7;AA=-NN/(48*DD);BB=MM/(864*DD*DD);RL=PolynomialRing(Kz,'lam');lam=RL.gen();fmod=lam**8+18*AA*lam**4+108*BB*lam**2-27*AA*AA
bf,bi,bI=rd(P/'lattice_bases.txt')
def mat(family):return matrix(Kz,[[evalrat(c,z)for c in row]for row in family])
BF=mat(bf);BI=mat(bi);II=mat(bI);assert all(M.rank()==8 for M in[BF,BI,II])
# Independent Hess/weak-Popov intersection of the two exported lattices.
M=II*BI.inverse();den=lcm([a.denominator()for a in M.list()]);PM=matrix(Rz,8,[(den*a).numerator()for a in M.list()]);U=identity_matrix(Rz,8)
def pivot(row):return max((int(a.degree()),j)for j,a in enumerate(row)if a)
steps=0
while True:
 piv=[pivot(row)for row in PM.rows()];pair=next(((i,j)for i in range(8)for j in range(i)if piv[i][1]==piv[j][1]),None)
 if pair is None:break
 i,j=pair
 if piv[i][0]<piv[j][0]:i,j=j,i
 col=piv[i][1];s=-PM[i,col].leading_coefficient()/PM[j,col].leading_coefficient()*Rz.gen()**(piv[i][0]-piv[j][0]);PM.add_multiple_of_row(i,j,s);U.add_multiple_of_row(i,j,s);steps+=1
assert PM==U*matrix(Rz,8,[(den*a).numerator()for a in M.list()])and U.det()in GF(5)and U.det()!=0
rowdegrees=[pivot(r)[0]for r in PM.rows()];dim=sum(max(0,den.degree()-d+1)for d in rowdegrees);assert dim==0
print('C242_GENUS17_INDEPENDENT_POPOV_DIMENSION_ZERO',rowdegrees,'den',den.degree(),'steps',steps,time.monotonic()-t0,flush=True)
if '--readonly' not in sys.argv:(P/'popov_stage.json').write_text(json.dumps(dict(dimension=0,rowdegrees=[int(x)for x in rowdegrees],denominator_degree=int(den.degree()),steps=steps),indent=2)+'\n')
# Maximality of the integral lattices.  Degree44 finite different and degree4 at infinity
# follow independently from genus17 and tame inertia 1^2*3^2 over z=infinity.
disc=fmod.discriminant();df=disc*BF.det()**2;di=disc*BI.det()**2;assert df.denominator().degree()==0 and df.numerator().degree()==44 and di.denominator().degree()-di.numerator().degree()==4
for name,B in [('finite',BF),('infinite',BI)]:
 for i,row in enumerate(B.rows()):
  p=RL(list(row));T=matrix(Kz,[[(p*lam**j%fmod)[k]for k in range(8)]for j in range(8)]);cp=T.charpoly()
  if name=='finite':assert all(c.denominator().degree()==0 for c in cp)
  else:assert all(not c or c.numerator().degree()<=c.denominator().degree()for c in cp)
  print('C242_GENUS17_INDEPENDENT_INTEGRAL_BASIS_ELEMENT',name,i+1,time.monotonic()-t0,flush=True)
print('C242_GENUS17_INDEPENDENT_MAXIMAL_BASES_CERTIFIED',time.monotonic()-t0,flush=True)
# Independent complete kernel11 isogenies, quotient coordinates, and actual norm support.
src,srcvals,deltavals=rd(P/'source_data.txt');K,E,P2,Q=curve(src);values=project17(E,P2,Q);delta=project17(E,P2,11*Q);assert values==[ev(K,a)for a in srcvals]and delta==[ev(K,a)for a in deltavals]and degree(values,9)==degree(delta,9)==3
assert project17(E,21*Q,P2+37*Q)==values;images=rd(P/'explicit_hecke_data.txt');RK=PolynomialRing(K,'w');normprod=RK.one();positive={};counts=Counter()
def key(vs):
 keys=[];vv=vs
 for j in range(vs[0].parent().degree()):keys.append(tuple(int(a)for v in vv for a in v._vector_()));vv=[v**5 for v in vv]
 return min(keys)
for i,(e,srcL,gdat,outL,embed,expect)in enumerate(images):
 L,EL,PL,QL=curve(srcL);LL,EG,PG,QG=curve(outL);assert L is LL and L.degree()==9 and e==1;zz=ev(L,embed);em=K.hom([zz],L);assert [em(a)for a in E.a_invariants()]==list(EL.a_invariants())and EL([em(a)for a in P2])==PL and EL([em(a)for a in Q])==QL
 RR=PolynomialRing(L,'w');g=RR([ev(L,a)for a in gdat]);iso=EL.isogeny(g,check=True);assert iso.degree()==11 and iso.codomain()==EG and iso(PL)==PG and iso(QL)==QG;vv=project17(EG,PG,QG);assert vv==[ev(L,a)for a in expect]and degree(vv,9)==9 and degree([vv[0]],9)==9
 B=matrix(GF(5),[(zz**j)._vector_()for j in range(9)]);cs=[]
 for a in g:
  co=B.solve_left(vector(GF(5),a._vector_()));b=ev(K,list(co));assert em(b)==a;cs.append(b)
 canon=[ev(K,list(B.solve_left(vector(GF(5),v._vector_()))))for v in vv];assert [em(v)for v in canon]==vv;kk=key(canon);positive[kk]=canon;counts[kk]+=1
 normprod*=RK(cs)
assert len(images)==12 and len(positive)==4 and set(counts.values())=={3}and normprod==RK(E.division_polynomial(11)).monic()
# The finite ideal is verified by containment and exact norm, using ordinary-place jets.
Tp=Rz.one()
for vs in positive.values():Tp*=Rz(vs[0].minpoly())
assert Tp.degree()==36 and Tp.is_squarefree();C=Tp*II*BF.inverse();assert all(a.denominator().degree()==0 for a in C.list());det=C.det();assert det.denominator().degree()==0 and det.numerator().degree()==288
polI=[RL(list(row))for row in II.rows()];polTI=[RL([Tp*c for c in row])for row in II.rows()]
def coeffeval(c,zz):return sum(zz.parent()(a)*zz**i for i,a in enumerate(c.numerator()))/sum(zz.parent()(a)*zz**i for i,a in enumerate(c.denominator()))
for vs in positive.values():
 zz,xx,yy=vs;uu=evalrat(udata,zz);tt=uu*uu;D0=tt*tt+13*tt+49;A0=-(tt*tt+5*tt+1)/(48*D0);ll=(3*xx*xx+A0)/(2*yy);R=PolynomialRing(zz.parent(),'l');ff=R([coeffeval(c,zz)for c in fmod]);assert ff.is_squarefree()and ff(ll)==0;other=ff.quo_rem(R([-ll,1]))[0]
 for b in polTI:assert R([coeffeval(c,zz)for c in b])%other==0
for vs,mult in[(values,1),(delta,11)]:
 zz,xx,yy=vs;L=zz.parent();PS=PowerSeriesRing(L,'e',default_prec=mult);ep=PS.gen();zs=PS(zz)+ep;uu=evalrat(udata,zs);tt=uu*uu;D0=tt*tt+13*tt+49;N0=tt*tt+5*tt+1;M0=tt**4+14*tt**3+63*tt*tt+70*tt-7;A0=-N0/(48*D0);B0=M0/(864*D0*D0);ac=A0[0];ll=PS((3*xx*xx+ac)/(2*yy))
 for j in range(5):
  f=ll**8+18*A0*ll**4+108*B0*ll*ll-27*A0*A0;fp=8*ll**7+72*A0*ll**3+216*B0*ll;assert fp[0]!=0;ll=(ll-f/fp).add_bigoh(mult)
 assert (ll**8+18*A0*ll**4+108*B0*ll*ll-27*A0*A0).valuation()>=mult
 for b in polI:
  val=sum(coeffeval(c,zs)*ll**j for j,c in enumerate(b));assert val.valuation()>=mult
 print('C242_GENUS17_INDEPENDENT_IDEAL_JET_PASS',mult,time.monotonic()-t0,flush=True)
out=dict(complete=True,rational_parameter_identity=True,prime=5,genus=17,source_orbit=8,quotient_closed_degree=3,independent_checked_11_isogenies=12,finite_integral_basis_discriminant_degree=44,infinite_integral_basis_discriminant_valuation=4,finite_ideal_norm_index=288,ordinary_zero_jet_orders=[1,11],independent_RR_dimension=0,popov_rowdegrees=[int(x)for x in rowdegrees],popov_denominator_degree=int(den.degree()))
if '--readonly'in sys.argv:assert json.loads((P/'independent_lattice_certificate.json').read_text())==out
else:(P/'independent_lattice_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C242_GENUS17_INDEPENDENT_LATTICE_SAGE_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
