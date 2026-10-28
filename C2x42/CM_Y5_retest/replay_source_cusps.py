from pathlib import Path
from math import gcd
from itertools import combinations_with_replacement
from collections import Counter
import json,sys
P=Path(__file__).resolve().parents[1]/'data'/'CM_Y5_retest';n=42;m=2
labels={(a,b,c,d)for a in range(n)for b in range(n)if gcd(gcd(a,b),n)==1 for c in range(2)for d in range(2)if(c*b-d*a)%2==1}
def orbit(x):
 a,b,c,d=x
 return frozenset(((s*(a+k*b))%n,s*b%n,(s*(c+k*d))%2,s*d%2)for s in[-1,1]for k in range(n))
classes=sorted(set(map(orbit,labels)),key=lambda c:sorted(c));idx={x:i for i,c in enumerate(classes)for x in c};assert len(classes)==72
def perm(fn):
 out=[]
 for c in classes:
  s={idx[fn(*x)]for x in c};assert len(s)==1;out.append(s.pop())
 assert len(set(out))==72;return out
gal={t:perm(lambda a,b,c,d:(t*a%n,b,t*c%2,d))for t in range(1,n)if gcd(t,n)==1};frob=gal[5]
def cycles(permutation):
 left=set(range(72));out=[]
 while left:
  start=min(left);cs=[];a=start
  while a not in cs:cs.append(a);left.remove(a);a=permutation[a]
  assert a==start;out.append(cs)
 return out
cc=cycles(frob);glob={frozenset(p[i]for p in gal.values())for i in range(72)};assert Counter(map(len,glob))=={1:18,2:9,6:6};assert Counter(map(len,cc))=={1:18,2:9,3:6,6:3}
A=perm(lambda a,b,c,d:((a+21*c)%42,(b+21*d)%42,a%2,b%2));B=perm(lambda a,b,c,d:(25*a%42,25*b%42,c,d));sub={'A3':A,'diamond7':B,'diag_plus':[A[B[i]]for i in range(72)],'diag_minus':[A[B[B[i]]]for i in range(72)]}
for cy in cc:
 globalorbit={p[cy[0]]for p in gal.values()}
 if len(cy)<=2:assert len(globalorbit)==len(cy)
 if len(cy)==3:
  assert len(globalorbit)==6 and {B[i]for i in cy}==set(cy)
  assert all(b%21==0 for i in cy for a,b,c,d in classes[i])
  for t,p in gal.items():
   ch=(1 if t%3==1 else-1)*(1 if t%7 in[1,2,4]else-1);assert({p[i]for i in cy}==set(cy))==(ch==1)
rat=[c[0]for c in cc if len(c)==1];quad=[c for c in cc if len(c)==2];cub=[c for c in cc if len(c)==3];divs=[Counter(a)for a in combinations_with_replacement(rat,3)]+[Counter([r]+q)for r in rat for q in quad]+[Counter(c)for c in cub];assert len(divs)==1308
stabs={h:sum(Counter({p[i]:v for i,v in D.items()})==D for D in divs)for h,p in sub.items()};assert stabs=={'A3':6,'diamond7':12,'diag_plus':6,'diag_minus':6}
out=dict(complete=True,geometric_cusps=72,global_signature={'1':18,'2':9,'6':6},F5_signature={'1':18,'2':9,'3':6,'6':3},all_degree_le_two_cusps_globally_intact=True,degree3_complements=1308,nonrational_degree3_complements=6,nonrational_cubic_cusps_b_stable=True,nonrational_cubic_field_quadratic_character=21,toric_7_and_3_cusp_levels=True,stable_cubic_complements=stabs)
if '--readonly'in sys.argv:assert json.loads((P/'source_cusp_certificate.json').read_text())==out
else:(P/'source_cusp_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C242_SOURCE_CUSP_GLOBAL_FROBENIUS_STABILIZER_REPLAY_PASS',json.dumps(out),flush=True)
