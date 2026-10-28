from pathlib import Path
import ast,json,sys
from collections import Counter
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
P=Path(__file__).resolve().parents[1]/'data'/'C2x40'
def data(name):return ast.literal_eval((P/name).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def fkey(pair,D):
 out=[]
 for j in range(D):out.append(tuple(int(x) for t in pair for x in t._vector_()));pair=[t**3 for t in pair]
 return min(out)
fp,cd,nc=data('source_frobenius.txt');assert len(fp)==115;tr=[];closed=[]
for n in range(1,10):
 s=-n*fp[114-n]-sum(fp[114-j]*tr[n-j-1] for j in range(1,n));tr.append(s);num=3**n+1-s-sum(d*closed[d-1] for d in range(1,n) if n%d==0);assert num%n==0;closed.append(num//n)
assert [v-Counter(cd)[d] for d,v in enumerate(closed,1)]==nc==[0,0,0,32,96,64,480,504,2048]
R=PolynomialRing(GF(3),'a');counts=[];orbitcounts=[]
for D in range(4,10):
 a=data('M2N40_D9_direct_full%d_arithmetic.txt'%D);assert a[0]==D and a[2];pol=a[1];rs=a[8];K=GF(3**D,'a',modulus=R(pol));z=K.gen();seen=set();countsby={int(d):0 for d in divisors(D)};own=0
 for i,(d,key) in enumerate(rs,1):
  b=sum(K(k)*z**j for j,k in enumerate(key[:D]));c=sum(K(k)*z**j for j,k in enumerate(key[D:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b)
  assert 2*T==E(0) and 40*Q==E(0) and 20*Q!=E(0) and 8*Q!=E(0) and T!=20*Q
  deg=min(int(e) for e in divisors(D) if b**(3**int(e))==b and c**(3**int(e))==c);assert deg==d;R2=20*Q;ims=set()
  for bb in range(2):
   TT=T+bb*R2
   for cc in range(2):
    for aa in [1,3,7,9,11,13,17,19]:
     U=cc*T+aa*Q;x=U[0]-TT[0];y=U[1];pair=[x**3/y**2,(c+3*TT[0])*x*x/y**2];ims.add(fkey(pair,D))
  assert tuple(key)==min(ims) and seen.isdisjoint(ims);seen.update(ims);countsby[d]+=len(ims);own+=d==D
 assert [countsby[int(d)] for d in divisors(D)]==a[7] and all(countsby[int(d)]==nc[int(d)-1] for d in divisors(D));counts.append(countsby[D]);orbitcounts.append(own)
assert orbitcounts==[1,3,2,15,21,64]
out=dict(complete=True,degree=9,prime=3,source_genus=57,closed_place_counts=nc,cusp_signature=dict(Counter(cd)),full_orbit_counts=orbitcounts,pure_orbits=sum(orbitcounts),full_action_order=32)
if '--readonly' in sys.argv:assert json.loads((P/'source_sage_certificate.json').read_text())==json.loads(json.dumps(out))
else:(P/'source_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C240_SAGE_COMPLETE_SOURCE_ENUMERATION_PASS',json.dumps(out),flush=True)
