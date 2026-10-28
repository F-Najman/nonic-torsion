from pathlib import Path
import ast,json,sys
from collections import Counter
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
P=Path(__file__).resolve().parent/'data'
def data(name):return ast.literal_eval((P/name).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def fkey(pair,D):
 out=[]
 for j in range(D):out.append(tuple(int(x) for t in pair for x in t._vector_()));pair=[t**5 for t in pair]
 return min(out)
fp,cd,nc=data('source_frobenius.txt');assert len(fp)==123;tr=[];closed=[]
for n in range(1,10):
 s=-n*fp[122-n]-sum(fp[122-j]*tr[n-j-1] for j in range(1,n));tr.append(s);num=5**n+1-s-sum(d*closed[d-1] for d in range(1,n) if n%d==0);assert num%n==0;closed.append(num//n)
assert [v-Counter(cd)[d] for d,v in enumerate(closed,1)]==nc==[0,0,0,126,648,2080,11376,48816,211776]
R=PolynomialRing(GF(5),'a');counts=[];orbitcounts=[]
for D in range(4,10):
 a=data('M2N42_D9_direct_full%d_arithmetic.txt'%D);assert a[0]==D and a[2];pol=a[1];rs=a[8];K=GF(5**D,'a',modulus=R(pol));z=K.gen();seen=set();countsby={int(d):0 for d in divisors(D)};own=0
 for i,(d,key) in enumerate(rs,1):
  b=sum(K(k)*z**j for j,k in enumerate(key[:D]));c=sum(K(k)*z**j for j,k in enumerate(key[D:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b)
  assert 2*T==E(0) and 42*Q==E(0) and 21*Q!=E(0) and 14*Q!=E(0) and 6*Q!=E(0) and T!=21*Q
  deg=min(int(e) for e in divisors(D) if b**(5**int(e))==b and c**(5**int(e))==c);assert deg==d;twos=[T,21*Q,T+21*Q];V=22*Q;ims=set()
  for ip in range(3):
   for iq in range(3):
    if ip==iq:continue
    for aa in [1,2,4,5,8,10]:
     U=twos[iq]+aa*V;TT=twos[ip];x=U[0]-TT[0];y=U[1];pair=[x**3/y**2,(c+3*TT[0])*x*x/y**2];ims.add(fkey(pair,D))
  assert tuple(key)==min(ims) and seen.isdisjoint(ims);seen.update(ims);countsby[d]+=len(ims);own+=d==D
 assert [countsby[int(d)] for d in divisors(D)]==a[7] and all(countsby[int(d)]==nc[int(d)-1] for d in divisors(D));counts.append(countsby[D]);orbitcounts.append(own)
assert len(orbitcounts)==6
out=dict(complete=True,degree=9,prime=5,source_genus=61,closed_place_counts=nc,cusp_signature=dict(Counter(cd)),full_orbit_counts=orbitcounts,pure_orbits=sum(orbitcounts),full_action_order=36)
if '--readonly' in sys.argv:assert json.loads((P/'source_sage_certificate.json').read_text())==json.loads(json.dumps(out))
else:(P/'source_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C242_SAGE_COMPLETE_SOURCE_ENUMERATION_PASS',json.dumps(out),flush=True)
