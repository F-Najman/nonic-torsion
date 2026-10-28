from pathlib import Path
import ast,json,re,sys,time
from collections import Counter
from sage.all import GF,QQ,ZZ,PolynomialRing,EllipticCurve,divisors,matrix,vector
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
S=Path(__file__).resolve().parents[1]/'data';P=S/'genus5_quotient';t0=time.monotonic()
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
# Rational descent identities, checked over Q(v), independent of both native finite probes.
RV=PolynomialRing(QQ,'v');V=RV.fraction_field();v=V.gen();b=v*v*(v-1);c=v*(v-1);H=v*v-v+1;s=(v**3-3*v+1)/(v*(v-1));t=s-8;D=t*t+13*t+49;N=t*t+5*t+1;M=t**4+14*t**3+63*t*t+70*t-7
ET=EllipticCurve(V,[1-c,-b,-b,0,0]);assert D==H**3/(v*v*(v-1)**2)and ET.c4()/H**4==N/D and ET.c6()/H**6==-M/(D*D)and ET.discriminant()/H**12==t/D**4
assert ET.discriminant()==t*v**8*(v-1)**8
assert ((v-1)/v)**2-((v-1)/v)+1==H/v**2
print('C242_G5_Q_DESCENT_IDENTITIES_SAGE_PASS',flush=True)

from sage.all import FunctionField
B0=FunctionField(QQ,'u');u=B0.gen();t=u*u;D=t*t+13*t+49;N=t*t+5*t+1;M=t**4+14*t**3+63*t*t+70*t-7;A=-N/(48*D);B=M/(864*D*D)
RX=PolynomialRing(B0,'x');F2=B0.extension(RX.gen()**3+A*RX.gen()+B,'r');r0=F2.gen();vr=F2(u/(4*D*D));r1=(-r0+vr/(3*r0*r0+F2(A)))/2;r2=-r0-r1;assert r0*r1+r0*r2+r1*r2==F2(A)and r0*r1*r2==-F2(B)
s=t+8;RV=PolynomialRing(F2,'v');v=RV.gen();F9=F2.extension(RV([F2(1),F2(s-3),-F2(s),F2(1)]),'v');v0=F9.gen();v1=v0*v0-F9(F2(s))*v0+F9(F2(s-2));v2=-v0*v0+F9(F2(s-1))*v0+2
w=F9(r0)*v0+F9(r1)*v1+F9(r2)*v2;wr=F9(r0)*v0+F9(r1)*v2+F9(r2)*v1
assert w**3-F9(N/48)*w+F9(M*(2*t+13)/(1728*D))-F9(u/(8*D))==0
assert wr**3-F9(N/48)*wr+F9(M*(2*t+13)/(1728*D))+F9(u/(8*D))==0
print('C242_DIAGONAL_RESOLVENT_Q_FUNCTIONFIELD_IDENTITY_PASS',flush=True)
raw=rd(S/'filter_quotient14_linear_arithmetic.txt');rows=[r for r in raw[2]if r[0]==9];fd=[a for a in raw[1]if a[0]==9][0];K=fld(fd[1]);checks=[]
for ii in [2,5,8]:
 row=rows[ii-1];b=ev(K,row[2][:9]);c=ev(K,row[2][9:]);E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P2=E(0,0);Q=E(b,b);S7=6*Q;x,y=S7[0],S7[1];m=(3*x*x+2*c*x+E.a4())/(2*y);a2p=3*x+c-m*m;a3p=2*y;kv=a3p/a2p;bv=-a2p**3/a3p**2;cv=1-2*m*a2p/a3p;v0=bv/cv;H=v0*v0-v0+1;r0=(P2[0]+c/3)/(kv*kv*H*H);r1=((21*Q)[0]+c/3)/(kv*kv*H*H);r2=((P2+21*Q)[0]+c/3)/(kv*kv*H*H);v1=(v0-1)/v0;v2=1/(1-v0);w=r0*v0+r1*v1+r2*v2;wr=r0*v0+r1*v2+r2*v1
 u,x3,y3=projection(E,P2,Q);t=u*u;D=t*t+13*t+49;N=t*t+5*t+1;M=t**4+14*t**3+63*t*t+70*t-7
 assert w**3-N/48*w+M*(2*t+13)/(1728*D)-u/(8*D)==0 and wr**3-N/48*wr+M*(2*t+13)/(1728*D)+u/(8*D)==0
 ds=[degree([u,x3,y3,w],9),degree([u,x3,y3,wr],9)];assert ds==([9,3]if ii==2 else[3,9]);checks.append(dict(source_orbit=ii,quotient_degrees=ds))
out=dict(complete=True,characteristic_zero_resolvent_identity=True,finite_source_checks=checks)
if '--readonly' in sys.argv:assert json.loads((P/'diagonal_resolvent_sage_certificate.json').read_text())==out
else:(P/'diagonal_resolvent_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C242_DIAGONAL_RESOLVENT_INDEPENDENT_SAGE_PASS',json.dumps(out),flush=True)
