from pathlib import Path
import ast,json,sys
from collections import Counter
from math import gcd,lcm
from sage.all import GF,PolynomialRing,EllipticCurve
P=Path(__file__).resolve().parent/'data'
def data(n):return ast.literal_eval((P/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def key(v,D):
 keys=[]
 for i in range(D):keys.append(tuple(int(x) for y in v for x in y._vector_()));v=[x**5 for x in v]
 return min(keys)
f=data('filter_quotient14_linear_arithmetic.txt');remain=data('filter_quotient14_isogeny_arithmetic.txt')[3];own={(r[0],r[1]):r for r in f[2]};sm={(r[0],r[1]):r for r in f[3]};R=PolynomialRing(GF(5),'z');fields={d:GF(5**d,'z',modulus=R(pol)) for d,pol in f[1]};aa=[1,2,4,5,8,10];autos=[]
for ip in range(3):
 for iq in range(3):
  if ip==iq:continue
  perm=[ip,iq,3-ip-iq];porder=None
  pp=list(range(3));porder=0
  while True:
   pp=[perm[i] for i in pp];porder+=1
   if pp==list(range(3)):break
  for a in aa:
   uorder=next(n for n in range(1,7) if pow(a,n,21) in [1,20]);order=lcm(porder,uorder)
   if order%2==0:autos.append((ip,iq,a,order))
assert len(autos)==27
curves={};cache={}
def im(row,h):
 d,i,ky=row[:3];kk=(d,tuple(ky))
 if kk not in curves:
  K=fields[d];z=K.gen();b=sum(K(k)*z**j for j,k in enumerate(ky[:d]));c=sum(K(k)*z**j for j,k in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b);curves[kk]=(E,c,[T,21*Q,T+21*Q],22*Q)
 ck=(kk,h)
 if ck in cache:return cache[ck]
 E,c,twos,V=curves[kk];ip,iq,a,o=h;TT=twos[ip];U=twos[iq]+a*V;x=U[0]-TT[0];y=U[1];ans=(d,key([x**3/y**2,(c+3*TT[0])*x*x/y**2],d));cache[ck]=ans;return ans
records=[];left=[]
for desc in remain:
 rows=[own[desc[1],desc[2]]]
 if desc[0]==2:rows.append(sm[desc[3],desc[4]])
 ncdegree=sum(r[0] for r in rows)
 if ncdegree!=8:left.append(desc);continue
 orig=Counter((r[0],tuple(r[2])) for r in rows);stabs=[h for h in autos if Counter(im(r,h) for r in rows)==orig];records.append(dict(descriptor=desc,even_order_stabilizers=stabs))
 if not stabs:left.append(desc)
out=dict(complete=True,source_action_order=36,even_order_elements=27,input_supports=len(remain),degree8_supports_checked=len(records),support_stabilizers=records,remaining=left,excluded_by_stabilizer=len(remain)-len(left))
if '--readonly' in sys.argv:assert json.loads((P/'stabilizer_filter_sage.json').read_text())==json.loads(json.dumps(out)) and data('after_even_stabilizer_descriptors.txt')==left
else:(P/'stabilizer_filter_sage.json').write_text(json.dumps(out,indent=2)+'\n');(P/'after_even_stabilizer_descriptors.txt').write_text(str(left)+'\n')
print('C242_EVEN_STABILIZER_SUPPORT_REPLAY_PASS',json.dumps({k:v for k,v in out.items() if k not in ['support_stabilizers','remaining']}),flush=True)
