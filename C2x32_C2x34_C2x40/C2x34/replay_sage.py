from pathlib import Path
import ast,json,sys
from collections import Counter
from sage.all import GF,PolynomialRing,EllipticCurve
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parents[1]/'data'/'C2x34'

def data(n):return ast.literal_eval((P/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def fkey(vs,D):
 out=[]
 for i in range(D):out.append(tuple(int(x) for v in vs for x in v._vector_()));vs=[v**3 for v in vs]
 return min(out)
class Quotient:
 def __init__(self,R,hs):
  self.n=len(R[0])-1;self.H=hermite_normal_form(Matrix([r[:-1] for r in R+hs]).T);self.index=abs(int(self.H.det()))
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0;w=list(v[:-1])
  for j in range(self.n-1,-1,-1):
   q=w[j]//int(self.H[j,j])
   for i in range(j+1):w[i]-=q*int(self.H[i,j])
  return tuple(w)
def tate(E,Q):
 x,y=Q[0],Q[1];a2,a4=E.a2(),E.a4();aa3=2*y;aa2=3*x+a2;aa4=3*x*x+2*x*a2+a4;aa1=2*aa4/aa3;aa2-= (aa4/aa3)**2
 b=-aa2**3/aa3**2;c=-(aa1*aa2-aa3)/aa3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1
 return [(s-r)/t,t/(s*s-s-r+1)]
def iso2(E,T,Q):
 K=E.base_field();r=T[0];A=E.a2()+3*r;B=3*r*r+2*E.a2()*r+E.a4();x=Q[0]-r;u=(x+A+B/x,Q[1]*(1-B/x**2));EP=EllipticCurve(K,[0,-2*A,0,A*A-4*B,0]);return EP,EP(u)
group,source,small,coverage,images,tests,remaining,fp,cdsource=data('four_quotients_cached_arithmetic.txt');invs,rels,hg,cd,cc=group
assert invs==[4,4,4,8760,1191360,0] and remaining==[] and len(source)==71 and len(tests)==239
GG=Quotient(rels,[]);QH=Quotient(rels,hg);assert QH.index==128 and all(QH.key(v)==(0,)*5 for v in cc)
tr=[];closed=[]
for n in range(1,10):
 z=-n*fp[98-n]-sum(fp[98-j]*tr[n-j-1] for j in range(1,n));tr.append(z);num=3**n+1-z-sum(d*closed[d-1] for d in range(1,n) if n%d==0);assert num%n==0;closed.append(num//n)
nc=[v-Counter(cdsource)[d] for d,v in enumerate(closed,1)];assert nc==[0,0,0,36,48,48,144,738,2160] and Counter(cdsource)=={1:24,8:3}
assert source==data('source_coverage.txt')[0] and small==data('source_coverage.txt')[1] and coverage==data('source_coverage.txt')[2]
R=PolynomialRing(GF(3),'z');fields={};src={};smallmap={};cov={r[0]:r for r in coverage};orbits=[2,1,2,3,18,45]
for d,pol,inputs,labels,os in coverage:
 assert len(inputs)==nc[d-1] and len(set(map(tuple,inputs)))==len(inputs) and len(os)==orbits[d-4]
 assert set(labels)==set(range(1,len(os)+1));assert {tuple(k) for ks in os for k in ks}==set(map(tuple,inputs));fields[d]=GF(3**d,'z',modulus=R(pol))
for row in source:
 d,i,pol,ky=row;K=fields[d];z=K.gen();b=sum(K(k)*z**j for j,k in enumerate(ky[:d]));c=sum(K(k)*z**j for j,k in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b);assert 34*Q==E(0) and 17*Q!=E(0) and 2*Q!=E(0) and T!=17*Q
 twos=[T,17*Q,T+17*Q];V=18*Q;ims=set()
 for ip in range(3):
  for iq in range(3):
   if ip==iq:continue
   for a in range(1,9):
    U=twos[iq]+a*V;PP=twos[ip];x=U[0]-PP[0];y=U[1];ims.add(fkey([x**3/y**2,(c+3*PP[0])*x*x/y**2],d))
 assert ims==set(map(tuple,cov[d][4][i-1]));assert ims=={tuple(k) for k,label in zip(cov[d][2],cov[d][3]) if label==i};assert tuple(ky) in ims;src[d,i]=row
for row in small:smallmap[row[0],row[1]]=row
cached={};actualisogenies=0
for key,e,w,cl in images:
 variant,d=key[:2];ky=key[2:];K=fields[d];z=K.gen();b=sum(K(k)*z**j for j,k in enumerate(ky[:d]));c=sum(K(k)*z**j for j,k in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b)
 if variant==1:EP,U=E,Q
 elif variant==2:EP,U=E,T+Q
 elif variant==3:EP,U=iso2(E,T,Q);actualisogenies+=1
 else:assert variant==4;EP,U=iso2(E,T+17*Q,Q);actualisogenies+=1
 assert 34*U==EP(0) and 17*U!=EP(0) and 2*U!=EP(0);bc=tate(EP,U);de=min(j for j in range(1,d+1) if d%j==0 and all(t**(3**j)==t for t in bc));assert de==e and d%e==0 and w==d//e
 cached[tuple(key)]=(e,w,cl)
expected=[[1,r[0],r[1]] for r in source]
for D,E in [[4,4],[4,5]]:
 for i in range(1,3):
  for j in range(1,nc[E-1]+1):expected.append([2,D,i,E,j])
assert len(expected)==239
mapvariants=0;normedsummands=0
for k,desc,cd,checks,passed in tests:
 assert desc==expected[k-1];rows=[src[desc[1],desc[2]]]
 if desc[0]==2:rows.append(smallmap[desc[3],desc[4]])
 assert cd==9-sum(r[0] for r in rows) and not passed;assert [q[0] for q in checks]==list(range(1,len(checks)+1))
 for v,cl,flag in checks:
  summed=[0]*len(invs)
  for d,i,pol,key in rows:
   e,w,vec=cached[tuple([v,d]+key)];summed=[x+y for x,y in zip(summed,vec)];normedsummands+=1
  assert GG.key(cl)==GG.key(summed) and bool(flag)==(QH.key(cl)==(0,)*5);mapvariants+=1
 assert not checks[-1][2] and all(q[2] for q in checks[:-1])
out=dict(complete=True,degree=9,excluded_groups=[[2,34],[1,68]],source_counts=nc,source_orbit_counts=orbits,full_action_order=48,source_places=sum(nc),supports=239,singles=71,mixed44=72,mixed45=96,checked_map_variants=mapvariants,normed_summands=normedsummands,distinct_actual_images=len(images),independent_velu2_images=actualisogenies,target_index=128,remaining=0)
if '--readonly' in sys.argv:assert json.loads((P/'independent_certificate.json').read_text())==out
else:(P/'independent_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C234_SAGE_COMPLETE_FINITE_SUPPORT_REPLAY_PASS',json.dumps(out),flush=True)
