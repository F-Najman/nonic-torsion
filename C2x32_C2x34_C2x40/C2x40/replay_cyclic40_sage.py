from pathlib import Path
import ast,json,sys,os,resource
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
resource.setrlimit(resource.RLIMIT_CPU,(570,600))
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parents[1]/'data'/'C2x40'
def data(n):return ast.literal_eval((P/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
class Quotient:
 def __init__(self,R,H):self.n=len(R[0])-1;self.M=hermite_normal_form(Matrix([r[:-1] for r in R+H]).T);self.index=abs(int(self.M.det()))
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0;w=list(v[:-1])
  for j in range(self.n-1,-1,-1):
   q=w[j]//int(self.M[j,j])
   for i in range(j+1):w[i]-=q*int(self.M[i,j])
  return tuple(w)
def tate(E,Q):
 x,y=Q[0],Q[1];aa3=2*y;aa2=3*x+E.a2();aa4=3*x*x+2*x*E.a2()+E.a4();aa1=2*aa4/aa3;aa2-=(aa4/aa3)**2;b=-aa2**3/aa3**2;c=-(aa1*aa2-aa3)/aa3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return [(s-r)/t,t/(s*s-s-r+1)]
f=data('filter_quotient20_linear_arithmetic.txt');own={(r[0],r[1]):r for r in f[2]};sm={(r[0],r[1]):r for r in f[3]};expected=data('filter_quotient20_isogeny_arithmetic.txt')[3];assert len(expected)==27
R=PolynomialRing(GF(3),'z');fields={d:GF(3**d,'z',modulus=R(pol)) for d,pol in f[1]};g,images,tests,remaining=data('filter_cyclic40_arithmetic.txt');iv,rels,hg,cd,cc,nativecomps=g
assert iv==[2,2,2,2,2,16,68880,895440,0] and remaining==[] and len(tests)==27;G=Quotient(rels,[]);H=Quotient(rels,hg);assert G.index==31579088486400 and H.index==4592
zero=(0,)*8;comps=[{zero},set(),set(),set()]
for d,v in zip(cd,cc):
 if d>3:continue
 for n in range(d,4):comps[n].update(H.key([a+b for a,b in zip(list(w)+[0],v)]) for w in list(comps[n-d]))
assert list(map(len,comps))==[1,1,5,5] and list(map(len,nativecomps))==[1,1,5,5]
cache={}
for kk,e,w,cl in images:
 variant,d=kk[:2];assert variant==1;ky=kk[2:];K=fields[d];z=K.gen();b=sum(K(a)*z**j for j,a in enumerate(ky[:d]));c=sum(K(a)*z**j for j,a in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P2=E(0,0);Q=E(b,b);assert 40*Q==E(0) and 20*Q!=E(0) and 8*Q!=E(0) and P2!=20*Q;xy=tate(E,Q);de=min(int(t) for t in divisors(d) if all(v**(3**int(t))==v for v in xy));assert de==e and w==d//e;cache[tuple(kk)]=cl
summands=0
for k,desc,cd,checks,passed in tests:
 assert desc==expected[k-1] and not passed and len(checks)==1 and checks[0][0]==1 and not checks[0][2];rows=[own[desc[1],desc[2]]]
 if desc[0]==2:rows.append(sm[desc[3],desc[4]])
 assert cd==9-sum(r[0] for r in rows) and cd in [0,1];v,cl,flag=checks[0];actual=[0]*len(iv)
 for row in rows:actual=[a+b for a,b in zip(actual,cache[tuple([1,row[0]]+row[2])])];summands+=1
 assert G.key(actual)==G.key(cl) and bool(flag)==(H.key([-x for x in cl]) in comps[cd])
out=dict(complete=True,degree=9,excluded_groups=[[2,40],[1,80]],target='X1(40)',target_genus=25,finite_Jacobian_order=G.index,rational_image_index=H.index,supports=27,cached_normalized_source_images=len(images),normed_summands=summands,checked_forgetful_maps=27,remaining=0,all_cases_fail_first_map=True,cusp_complement_counts=[1,1,5,5])
if '--readonly' in sys.argv:assert json.loads((P/'cyclic40_sage_certificate.json').read_text())==out
else:(P/'cyclic40_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C240_CYCLIC40_INDEPENDENT_FINITE_REPLAY_PASS',json.dumps(out),flush=True)
