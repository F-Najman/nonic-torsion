from pathlib import Path
import ast,json,sys
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent/'data'

def data(n):return ast.literal_eval((P/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
class Quotient:
 def __init__(self,R,hs):
  self.n=len(R[0])-1;self.H=hermite_normal_form(Matrix([r[:-1] for r in R+hs]).T);self.index=abs(int(self.H.det()))
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0;w=list(v[:-1])
  for j in range(self.n-1,-1,-1):
   q=w[j]//int(self.H[j,j])
   for i in range(j+1):w[i]-=q*int(self.H[i,j])
  return tuple(w)
g,fields,single,small,pairs,triples,survivors,coverage=data('filter_quotient18_linear_arithmetic.txt');invs,rels,cd,cc,hg,hinv=g
assert invs==[6,84,252,0] and hinv==[2,42,126];QH=Quotient(rels,hg);GG=Quotient(rels,[]);assert QH.index==12 and all(QH.key(v)==(0,0,0) for v in cc)
assert len(single)==11800 and len(small)==3792 and len(pairs)==11784 and len(triples)==2352 and len(survivors)==2064
R=PolynomialRing(GF(5),'z');ks={d:GF(5**d,'z',modulus=R(pol)) for d,pol in fields};own={};sm={};degreemaps=0
for is_small,rows in [(False,single),(True,small)]:
 for row in rows:
  d,i,key,e,w,cl=row[:6];K=ks[d];z=K.gen();b=sum(K(t)*z**j for j,t in enumerate(key[:d]));c=sum(K(t)*z**j for j,t in enumerate(key[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);Q=E(b,b);U=2*Q;x,y=U[0],U[1];bc=[x**3/y**2,c*x*x/y**2];de=min(int(t) for t in divisors(d) if all(v**(5**int(t))==v for v in bc));assert e==de and w==d//e
  (sm if is_small else own)[d,i]=row;degreemaps+=1
  if not is_small:assert bool(row[6])==(QH.key(cl)==(0,0,0))
expected=[]
for d,i in own:
 if own[d,i][6]:expected.append([1,d,i])
for d,es in [(3,3),(3,4),(3,5),(3,6),(4,4),(4,5)]:
 rs=[r for r in pairs if r[:1]==[d] and r[2]==es];assert len(rs)==sum(k[0]==d for k in own)*sum(k[0]==es for k in sm)
 assert [(r[1],r[3]) for r in rs]==[(i,j) for dd,i in own if dd==d for ee,j in sm if ee==es]
 for D,i,E,j,cl,yes in rs:
  actual=[a+b for a,b in zip(own[D,i][5],sm[E,j][5])];assert GG.key(actual)==GG.key(cl) and bool(yes)==(QH.key(cl)==(0,0,0))
  if yes:expected.append([2,D,i,E,j])
assert [(r[1],r[2],r[3]) for r in triples]==[(i,j,k) for d,i in own if d==3 for e,j in sm if e==3 for f,k in sm if f==3 and k>=j]
for D,i,j,k,cl,yes in triples:
 actual=[a+b+c for a,b,c in zip(own[3,i][5],sm[3,j][5],sm[3,k][5])];assert GG.key(actual)==GG.key(cl) and bool(yes)==(QH.key(cl)==(0,0,0))
 if yes:expected.append([3,3,i,3,j,3,k])
assert expected==survivors
out=dict(complete=True,prime=5,target='X1(2,18)',target_image_index=12,single_supports=11800,pair_supports=11784,triple_supports=2352,total_supports=25936,normalized_source_image_degrees_checked=degreemaps,remaining=2064)
if '--readonly' in sys.argv:assert json.loads((P/'first_filter_sage_certificate.json').read_text())==out
else:(P/'first_filter_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C236_SAGE_COMPLETE_FIRST_FILTER_REPLAY_PASS',json.dumps(out),flush=True)
