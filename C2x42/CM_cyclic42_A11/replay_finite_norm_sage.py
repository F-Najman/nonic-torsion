from pathlib import Path
import ast,json,re,sys,time
from collections import Counter
from sage.all import GF,QQ,ZZ,PolynomialRing,EllipticCurve,divisors,matrix,vector
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parents[1]/'data'/'CM_cyclic42_A11';S=P.parent;ROOT=S.parent.parent;t0=time.monotonic()
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
proof=rd(P/'explicit_hecke_data.txt');assert [r[0]for r in proof]==[1,2,6,8];logs=(ROOT/'C2x42/logs/CM_cyclic42_A11/test.log').read_text();assert all(re.search(r'C242_CM_CYCLIC42_A11_CASE '+str(i)+r' 0 ',logs)for i in[1,2,6,8])
cm=rd(S/'CM_Y5_retest/retests.txt');assert cm[1]==[1,2,6,8];assert all(r[4]==[0]for ii in[1,2,6,8]for r in cm[0][ii-1][3]);count=0;outputs=[]
for ii,src,images in proof:
 K,E,P2,Q=curve(src);assert K.degree()==6 and E.j_invariant()==0 and Q.order()==42 and P2.order()==2 and P2!=21*Q;orig=xy(E,Q);delta=xy(E,11*Q);d=degree(orig,6);assert degree(delta,6)==d and eq(*orig)==eq(*delta)==0;RK=PolynomialRing(K,'w');normprod=RK.one();normdeg=0
 for e,srcL,gdat,outL,embed,expect,de,w in images:
  L,EL,PL,QL=curve(srcL);LL,EG,PG,QG=curve(outL);assert L is LL and L.degree()==6 and e==1;zz=ev(L,embed);em=K.hom([zz],L);assert [em(a)for a in E.a_invariants()]==list(EL.a_invariants())and EL([em(a)for a in P2])==PL and EL([em(a)for a in Q])==QL
  RL=PolynomialRing(L,'w');g=RL([ev(L,a)for a in gdat]);iso=EL.isogeny(g,check=True);assert iso.degree()==11 and iso.codomain()==EG and iso(PL)==PG and iso(QL)==QG;vals=xy(EG,QG);assert vals==[ev(L,a)for a in expect]and eq(*vals)==0;assert degree(vals,6)==de and w==6//de;normdeg+=w*de
  B=matrix(GF(5),[(zz**j)._vector_()for j in range(6)]);co=[]
  for a in g:
   vec=B.solve_left(vector(GF(5),a._vector_()));b=ev(K,list(vec));assert em(b)==a;co.append(b)
  normprod*=RK(co);count+=1
 assert len(images)==12 and normprod==RK(E.division_polynomial(11)).monic()and normdeg==72;outputs.append(dict(source_orbit=ii,target_degree=d,norm_weight=6//d,Hecke_degree=72,native_RR_dimension=0))
out=dict(complete=True,prime=5,target='X1(42)',cases=outputs,independent_checked_11_isogenies=count,all_norm_multiplicities_verified=True,all_plane_equations_verified=True,globally_rational_complement_condition_verified=True,native_actual_exit_code=0)
if '--readonly'in sys.argv:assert json.loads((P/'finite_norm_certificate.json').read_text())==out
else:(P/'finite_norm_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C242_FOUR_CM_CYCLIC42_INDEPENDENT_FINITE_NORM_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
