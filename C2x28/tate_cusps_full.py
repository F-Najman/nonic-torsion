"""Independent exact Tate cusp permutations and rational cusp lifting at3."""
from pathlib import Path
from math import gcd
from collections import Counter
import json
P=Path(__file__).resolve().parent;n=28
labels={(a,b,c,d) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1 for c in range(2) for d in range(2) if (c*b-d*a)%2==1}
def orbit(x):
 a,b,c,d=x
 return frozenset(((s*(a+k*b))%n,s*b%n,s*(c+k*d)%2,s*d%2) for s in (-1,1) for k in range(n))
classes=sorted(set(map(orbit,labels)),key=lambda s:sorted(s));index={x:i for i,cl in enumerate(classes) for x in cl}
def galperm(t):
 ans=[]
 for cl in classes:
  images={index[(t*a%28,b,t*c%2,d)] for a,b,c,d in cl};assert len(images)==1;ans.append(images.pop())
 assert len(set(ans))==len(classes);return ans
perms={t:galperm(t) for t in range(28) if gcd(t,28)==1};frob=perms[3]
left=set(range(len(classes)));cycles=[]
while left:
 i=min(left);cy=[]
 while i not in cy:cy.append(i);left.remove(i);i=frob[i]
 assert i==cy[0];cycles.append(cy)
assert Counter(map(len,cycles))=={1:12,2:6,3:4,6:2}
low=[cy for cy in cycles if len(cy)<=2]
assert all(all(pp[cy[0]] in cy for pp in perms.values()) for cy in low)
# h(P,Q)=(P+14Q,13Q).
h=[]
for cl in classes:
 images={index[(13*a%28,13*b%28,(c+a)%2,(d+b)%2)] for a,b,c,d in cl};assert len(images)==1;h.append(images.pop())
assert all(h[h[i]]==i for i in range(len(h)))
rat=[cy[0] for cy in cycles if len(cy)==1];assert all(h[i]!=i for i in rat)
a=dict(level=28,prime=3,geometric_cusps=len(classes),closed_cusp_degrees=sorted(map(len,cycles)),rational_and_quadratic_residue_cusps_are_global_closed_cusps=True,stabilizer=[1,0,13],rational_cusps=rat,involution_permutation=h,fixed_rational_cusps=[],complete=True)
# stdout only: safe for read-only certificate replay.
print(json.dumps(a,indent=2))
print('C228_D9_TATE_CUSP_AND_CORRECTED_STABILIZER_PASS rational=12 quadratic=6 fixed_rational_cusps=0')
