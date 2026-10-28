#!/usr/bin/env python3
"""Independent integer-lattice replay of the completed F3 class packet."""
from pathlib import Path
from collections import Counter
from math import prod,comb
import json
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parents[2]/'C35/data/presentation_v2'
D=json.loads((P/'class_data_fb.json').read_text());assert D['complete']
R=sp.Matrix(D['relations']);dim=R.cols;free=[j for j in range(dim) if all(R[i,j]==0 for i in range(R.rows))]
assert len(free)==1;free=free[0];ix=[j for j in range(dim) if j!=free]
assert dim==6 and all(R[i,j]==0 for i in range(R.rows) for j in range(dim) if i!=j)
invs=[int(R[i,i]) for i in ix];assert invs==[4,4,8,520,218297040]
order=prod(invs);assert order==14529850982400
hq=D['rational_cusp_generators'];assert all(len(a)==dim and a[free]==0 for a in hq)
lat=sp.Matrix([[int(a[j]) for j in ix] for a in D['relations']+hq]);B=hermite_normal_form(lat.T)
assert B.shape==(5,5) and abs(B.det())==10240
BI=B.inv();assert order//10240==13*109148520

def canonical(a):
 assert len(a)==dim and a[free]==0
 return tuple(int(a[j])%n for j,n in zip(ix,invs))

def member(a):
 v=sp.Matrix([int(a[j]) for j in ix]);w=BI*v
 return all(x.q==1 for x in w)

known={canonical(a) for a in D['known_octic_classes']};assert len(known)==6 and all(member(a) for a in D['known_octic_classes'])
expected=[0,0,12,0,36,126,348,813,2240];assert D['noncuspidal_counts']==expected
f=D['frobenius_coefficients'];assert len(f)==51 and f[0]==3**25 and f[-1]==1 and sum(f)==order
assert all(f[i]==3**(25-i)*f[50-i] for i in range(26))
traces=[];closed=[]
for e in range(1,10):
 t=-sum(f[50-i]*traces[e-i-1] for i in range(1,e))-e*f[50-e];traces.append(t)
 cusps=sum(d*n for d,n in [(1,12),(4,3),(6,2),(12,1)] if e%d==0)
 nc=3**e+1-t-cusps
 v=nc-sum(j*closed[j-1] for j in range(1,e) if e%j==0)
 assert v%e==0;closed.append(v//e)
assert closed==expected
pure={};seen=set();types=Counter();unexpected=[];seen8=set()
for t in D['tests']:
 typ=tuple(t['degrees']);ids=tuple(t['indices']);key=typ,ids
 assert key not in seen;seen.add(key);types[typ]+=1
 assert len(typ)==len(ids) and sum(typ)<=9
 if len(typ)==1:pure[typ[0],ids[0]]=t['class']
 truth=member(t['class']);assert truth==t['in_H']
 is_known=sum(typ)==8 and canonical(t['class']) in known;assert is_known==t['known_octic']
 if truth:
  if is_known:seen8.add(canonical(t['class']))
  else:unexpected.append({'degrees':typ,'indices':ids,'class':t['class']})
for t in D['tests']:
 if len(t['degrees'])<2:continue
 a=[0]*dim
 for d,i in zip(t['degrees'],t['indices']):a=[x+y for x,y in zip(a,pure[d,i])]
 assert canonical(a)==canonical(t['class'])
assert types==Counter({(3,):12,(5,):36,(6,):126,(7,):348,(8,):813,(9,):2240,(3,3):78,(5,3):432,(6,3):1512,(3,3,3):364})
assert len(D['tests'])==5961
out={'order':order,'quotient_index':10240,'tests':5961,'frobenius_counts_replayed':True,'unexpected':unexpected,'known_octic_classes':len(seen8),'supports_exclusion':not unexpected and seen8==known}
print('NONIC35_INDEPENDENT_CLASS_REPLAY '+json.dumps(out))
