from pathlib import Path
from math import gcd
from collections import Counter
import json,sys
P=Path(__file__).resolve().parents[1]/'data'/'quotient9';N=57
U=[u for u in range(N)if gcd(u,N)==1];H=[];h=1
while h not in H:H.append(h);h=h*4%N
assert len(H)==9
def key(a,c):
 d=gcd(c,N);a%=d
 return min((a,c%N),((-a)%d,(-c)%N))
base={key(a,c)for c in range(N)for a in range(gcd(c,N))if gcd(a,gcd(c,N))==1}
def diamond(P,u):a,c=P;return key(a*pow(u,-1,N),c*u)
def orbH(P):return frozenset(diamond(P,u)for u in H)
cs={orbH(P)for P in base};assert len(cs)==8
index={v:o for o in cs for v in o}
def gal(o,u):a,c=min(o);return index[key(a*u,c)]
def orbits(U):
 todo=set(cs);ans=[]
 while todo:
  p=min(todo,key=lambda o:sorted(o));o=frozenset(gal(p,u)for u in U);assert o<=todo;todo-=o;ans.append(o)
 return ans
F=[];v=1
while v not in F:F.append(v);v=v*2%N
G=orbits(U);Q=orbits(F)
def width(o):return N//gcd(min(o)[1],N)
rows=[]
for w in [1,3,19,57]:
 rows.append(dict(width=w,global_degrees=sorted(len(o)for o in G if width(next(iter(o)))==w),finite_degrees=sorted(len(o)for o in Q if width(next(iter(o)))==w)))
out=dict(complete=True,level=N,H_generator=4,H_order=9,geometric_cusps=8,global_signature=dict(Counter(map(len,G))),F2_signature=dict(Counter(map(len,Q))),width_rows=rows)
print('NONIC57_GLOBAL_CUSP_AUDIT_PASS',json.dumps(out))
# The rational-q-expansion model is the dual moduli model. Fricke changes
# the displayed j-width w to N/w while preserving closed cusp fields.
out['dual_width_rows']=[dict(width=N//r['width'],global_degrees=r['global_degrees'],finite_degrees=r['finite_degrees'])for r in rows]
source=list(base)
def sourcegal(p,u):return key(p[0]*u,p[1])
def sourceorbits(U):
 todo=set(source);ans=[]
 while todo:
  p=min(todo);o=frozenset(sourcegal(p,u)for u in U);assert o<=todo;todo-=o;ans.append(o)
 return ans
SG=sourceorbits(U);SF=sourceorbits(F)
assert Counter(map(len,SG))=={1:18,2:9,18:2}
assert Counter(map(len,SF))=={1:18,2:9,9:2,18:1}
assert all(o in SG for o in SF if len(o)<=3)
bad=[o for o in SF if len(o)<=9 and o not in SG];assert len(bad)==2 and all(len(o)==9 and N//gcd(min(o)[1],N)==1 for o in bad)
out.update(source_global_signature=dict(Counter(map(len,SG))),source_F2_signature=dict(Counter(map(len,SF))),all_source_cusps_degree_le3_globally_intact=True,only_non_global_source_cusps_degree_le9_two_single_nonics=True,non_global_nonic_dual_target_width=57,non_global_nonic_norm_weight=9)
if '--readonly' in sys.argv:assert json.loads((P/'cusp57.json').read_text())==json.loads(json.dumps(out))
else:(P/'cusp57.json').write_text(json.dumps(out,indent=2)+'\n')
print('NONIC57_SOURCE_AND_DUAL_CUSP_AUDIT_PASS',json.dumps(out))
