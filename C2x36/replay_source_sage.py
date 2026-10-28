from pathlib import Path
import ast,json,sys,resource
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
resource.setrlimit(resource.RLIMIT_CPU,(14400,14430))
P=Path(__file__).resolve().parent/'data'

def data(name):return ast.literal_eval((P/name).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def fkey(pair,D):
 out=[]
 for j in range(D):out.append(tuple(int(x) for t in pair for x in t._vector_()));pair=[t**5 for t in pair]
 return min(out)
fp,nc=data('source_frobenius.txt');assert len(fp)==83
tr=[];closed=[]
for n in range(1,10):
 s=-n*fp[82-n]-sum(fp[82-j]*tr[n-j-1] for j in range(1,n));tr.append(s)
 num=5**n+1-s-sum(d*closed[d-1] for d in range(1,n) if n%d==0);assert num%n==0;closed.append(num//n)
assert [v-({1:24,2:8,6:4}.get(d,0)) for d,v in enumerate(closed,1)]==nc==[0,0,48,120,720,2904,10416,50334,218080]
R=PolynomialRing(GF(5),'a');counts=[];orbitcounts=[]
for D in range(3,10):
 a=data('M2N36_D9_direct_full%d_arithmetic.txt'%D);assert a[0]==D and a[2];pol=a[1];rs=a[8]
 K=GF(5**D,'a',modulus=R(pol));z=K.gen();seen=set();countsby={int(d):0 for d in divisors(D)};own=0
 for i,(d,key) in enumerate(rs,1):
  b=sum(K(k)*z**j for j,k in enumerate(key[:D]));c=sum(K(k)*z**j for j,k in enumerate(key[D:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b)
  assert 2*T==E(0) and 36*Q==E(0) and 18*Q!=E(0) and 12*Q!=E(0) and T!=18*Q
  deg=min(int(e) for e in divisors(D) if b**(5**int(e))==b and c**(5**int(e))==c);assert deg==d
  R2=18*Q;ims=set()
  for bb in range(2):
   TT=T+bb*R2
   for cc in range(2):
    for aa in [1,5,7,11,13,17]:
     U=cc*T+aa*Q;x=U[0]-TT[0];y=U[1];pair=[x**3/y**2,(c+3*TT[0])*x**2/y**2];ims.add(fkey(pair,D))
  assert tuple(key)==min(ims) and seen.isdisjoint(ims);seen.update(ims);countsby[d]+=len(ims);own+=d==D
  if i%500==0:print('C236_SAGE_SOURCE_PROGRESS',D,i,flush=True)
 assert [countsby[int(d)] for d in divisors(D)]==a[7]
 assert all(countsby[int(d)]==nc[int(d)-1] for d in divisors(D));counts.append(countsby[D]);orbitcounts.append(own)
 print('C236_SAGE_COMPLETE_SOURCE_DEGREE_PASS',D,own,countsby,flush=True)
assert orbitcounts==[2,5,30,126,434,2115,9088]
out=dict(complete=True,degree=9,prime=5,closed_place_counts=nc,full_orbit_counts=orbitcounts,pure_orbits=sum(orbitcounts),full_action_order=24)
if '--readonly' in sys.argv:assert json.loads((P/'source_sage_certificate.json').read_text())==out
else:(P/'source_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C236_SAGE_COMPLETE_SOURCE_ENUMERATION_PASS',flush=True)
