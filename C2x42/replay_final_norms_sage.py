from pathlib import Path
import ast,json,re,sys,time
from collections import Counter
from sage.all import GF,QQ,ZZ,PolynomialRing,EllipticCurve,divisors,matrix,vector
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent/'data';S=P.parent;ROOT=S.parent;t0=time.monotonic()
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

from sage.all import sage_eval
import re
Rxy=PolynomialRing(GF(5),2,names=('x','y'));xx,yy=Rxy.gens();eq=sage_eval((ROOT/'mdmagma/models_X1_n/FFFc42.txt').read_text().replace('\n',' '),locals={'x':xx,'y':yy})
def xy(E,Q):
 x,y=Q[0],Q[1];m=(3*x*x+2*E.a2()*x+E.a4())/(2*y);a2=3*x+E.a2()-m*m;a3=2*y;b=-a2**3/a3**2;c=1-2*m*a2/a3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return[(s-r)/t,t/(s*s-s-r+1)]

input=rd(P/'filter_quotient14_linear_arithmetic.txt');own={(r[0],r[1]):r for r in input[2]};small={(r[0],r[1]):r for r in input[3]};fieldpolys={d:pol for d,pol in input[1]};sourceproof=rd(P/'cyclic42_A11_seven/explicit_hecke_data.txt');altproof=rd(P/'quartic_alt_A11/explicit_hecke_data.txt');cases=[[1,4,5],[1,6,39],[2,4,3,4,55],[2,4,3,4,65],[2,4,3,5,260],[2,4,3,5,461],[2,4,3,5,589]]
def pk(vals):
 L=vals[0].parent();d=degree(vals,L.degree());a=None
 for aa in range(5):
  g=vals[0]+aa*vals[1]
  if degree([g],L.degree())==d:a=aa;break
 if a is None:g=vals[1];assert degree([g],L.degree())==d;a=5
 B=matrix(GF(5),[(g**j)._vector_()for j in range(d)]);cx=B.solve_left(vector(GF(5),vals[0]._vector_()));cy=B.solve_left(vector(GF(5),vals[1]._vector_()));pol=list(map(int,g.minpoly()));return(tuple(pol),tuple(map(int,cx)),tuple(map(int,cy)),a)
points={}
def geom(vals):
 k=pk(vals);d=len(k[0])-1;point=dict(polynomial=list(k[0]),x_coefficients=list(k[1]),y_coefficients=list(k[2]),primitive_linear_combination=k[3],degree=d,unramified_y=bool(eq.derivative(yy)(*vals)!=0),x_primitive=degree([vals[0]],vals[0].parent().degree())==d)
 sk=json.dumps(k,separators=(',',':'));points[sk]=point;return sk,d
count=0;weighted=0;Apoints={};sourcekeys=[]
def evaluate(src,images,expectedD):
 global count,weighted
 K,E,Paux,Q=curve(src);d=K.degree();assert d==expectedD and Q.order()==42;orig=xy(E,Q);delta=xy(E,11*Q);assert eq(*orig)==eq(*delta)==0;ko,de=geom(orig);kd,dde=geom(delta);assert de==dde and d%de==0;D=Counter();D[ko]-=d//de;D[kd]-=11*d//de;RK=PolynomialRing(K,'w');normprod=RK.one();normdegree=0
 for e,srcL,gdat,outL,embed,expect,dd,w in images:
  L,EL,PL,QL=curve(srcL);LL,EG,PG,QG=curve(outL);assert L is LL and L.degree()==d*e;zz=ev(L,embed);em=K.hom([zz],L);assert [em(a)for a in E.a_invariants()]==list(EL.a_invariants())and EL([em(a)for a in Q])==QL and EL([em(a)for a in Paux])==PL
  RL=PolynomialRing(L,'w');g=RL([ev(L,a)for a in gdat]);iso=EL.isogeny(g,check=True);assert iso.degree()==11 and iso.codomain()==EG and iso(PL)==PG and iso(QL)==QG;vals=xy(EG,QG);assert vals==[ev(L,a)for a in expect]and eq(*vals)==0;k,de=geom(vals);assert de==dd and d*e%de==0 and w==d*e//de;D[k]+=w;normdegree+=w*de
  ng=RL.one();gg=g
  for j in range(e):ng*=gg;gg=RL([a**(5**d)for a in gg])
  assert gg==g;B=matrix(GF(5),[(zz**j)._vector_()for j in range(d)]);co=[]
  for a in ng:
   vec=B.solve_left(vector(GF(5),a._vector_()));b=ev(K,list(vec));assert em(b)==a;co.append(b)
  normprod*=RK(co);count+=1;weighted+=e
 assert sum(r[0]for r in images)==12 and normprod==RK(E.division_polynomial(11)).monic()and normdegree==12*d
 D={k:int(v)for k,v in D.items()if v};assert sum(points[k]['degree']*v for k,v in D.items())==0;return D
for kk,src,de,w,images in sourceproof:
 d=kk[0];K,E,P2,Q=curve(src);assert src[0]==fieldpolys[d]and P2==E(0,0)and P2.order()==2 and P2!=21*Q;assert list(map(int,Q[0]._vector_()))==kk[1:1+d]and list(map(int,E.a2()._vector_()))==kk[1+d:]
 assert degree(xy(E,Q),d)==de and w==d//de;Apoints[tuple(kk)]=evaluate(src,images,d);sourcekeys.append(tuple(kk));print('C242_FINAL_SOURCE_SUMMAND_NORM_PASS',d,len(sourcekeys),flush=True)
expected=set();casegeoms={}
for j,desc in enumerate(cases,1):
 rows=[own[desc[1],desc[2]]]
 if desc[0]==2:rows.append(small[desc[3],desc[4]])
 D=Counter()
 for r in rows:
  kk=tuple([r[0]]+r[2]);expected.add(kk)
  for k,v in Apoints[kk].items():D[k]+=v
 casegeoms['seven_'+str(j)]={k:int(v)for k,v in D.items()if v}
assert expected==set(sourcekeys)and len(expected)==8
for mapno,src,target,images in altproof:
 K,E,P2,Q=curve(src);L,EG,Qa,Qb=curve(target);assert K is L and Qa==Qb
 if mapno==2:assert EG==E and Qa==P2+Q
 else:
  TT=P2 if mapno==3 else P2+21*Q;iso2=E.isogeny(TT,check=True);assert iso2.codomain()==EG and iso2(Q)==Qa
 casegeoms['alt_'+str(mapno)]=evaluate(target,images,4)
for ii,src,images in rd(P/'CM_cyclic42_A11/explicit_hecke_data.txt'):casegeoms['CM_'+str(ii)]=evaluate(src,images,6)
assert set(casegeoms)=={'seven_'+str(j)for j in range(1,8)}|{'alt_2','alt_3','alt_4','CM_1','CM_2','CM_6','CM_8'}
for folder,marker in[('cyclic42_A11_seven','C242_A11_SEVEN_DONE'),('quartic_alt_A11','C242_QUARTIC_ALT_A11_DONE'),('CM_cyclic42_A11','C242_CM_CYCLIC42_A11_DONE')]:
 assert marker in(S/'logs'/folder/'test.log').read_text()
for j in range(1,8):assert 'C242_A11_SEVEN_CASE_DONE '+str(j)+' '+str(1 if j==1 else 0)+' 'in(S/'logs/cyclic42_A11_seven/test.log').read_text()
assert all('C242_QUARTIC_ALT_A11_CASE_DONE '+str(j)+' 0 'in(S/'logs/quartic_alt_A11/test.log').read_text()for j in[2,3,4])
geomout=dict(points=points,cases=casegeoms,complete=True);out=dict(complete=True,prime=5,target='X1(42)',independent_kernel_orbits=count,geometric_kernel_count=weighted,source_summands=8,nonCM_cases=7,alternative_maps=3,CM_cases=4,all_plane_equations_and_norms=True,all_14_case_divisors_reconstructed=True,all_source_case_links=True,x_nonprimitive_points=sum(not p['x_primitive']for p in points.values()),ramified_y_points=sum(not p['unramified_y']for p in points.values()))
if '--readonly'in sys.argv:assert json.loads((P/'final_norm_certificate.json').read_text())==out and json.loads((P/'final_divisor_geometry.json').read_text())==geomout
else:(P/'final_norm_certificate.json').write_text(json.dumps(out,indent=2)+'\n');(P/'final_divisor_geometry.json').write_text(json.dumps(geomout,indent=2)+'\n')
print('C242_FINAL_DIRECT_HECKE_INDEPENDENT_FINITE_NORM_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
