#!/usr/bin/env python3
from collections import Counter
from math import gcd
from pathlib import Path
import json
P=Path(__file__).resolve().parent;n=133
labels={(a,b) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1}
def orbit(x):
 a,b=x
 return frozenset(((s*(a+k*b))%n,s*b%n) for s in (-1,1) for k in range(n))
classes=sorted(set(map(orbit,labels)),key=lambda c:sorted(c));assert len(classes)==216
index={x:i for i,c in enumerate(classes) for x in c};reps=[min(c) for c in classes];units=[u for u in range(n) if gcd(u,n)==1]
gal={u:[index[(u*a%n,b)] for a,b in reps] for u in units}
hs={pow(78,j,n) for j in range(18)};assert len(hs)==18 and all(h%7==1 for h in hs)
qos=sorted({tuple(sorted({index[(h*a%n,h*b%n)] for h in hs})) for a,b in reps});assert len(qos)==12 and all(len(o)==18 for o in qos)
qi={i:j for j,o in enumerate(qos) for i in o}
qgal={u:[qi[gal[u][o[0]]] for o in qos] for u in units}
assert all(qi[gal[u][i]]==qgal[u][j] for u in units for j,o in enumerate(qos) for i in o)
left=set(range(12));go=[]
while left:
 i=min(left);o={qgal[u][i] for u in units};left-=o;go.append(o)
left=set(range(12));fo=[]
while left:
 i=min(left);o=[]
 while i not in o:o.append(i);left.remove(i);i=qgal[2][i]
 assert i==o[0];fo.append(o)
assert Counter(map(len,go))==Counter(map(len,fo))=={1:6,3:2}
out={'level':133,'H_generator':78,'H_order':18,'H_fixes_P7':True,'source_cusp_count':216,'H_free_on_source_cusps':True,'quotient_cusp_count':12,'global_signature':dict(Counter(map(len,go))),'F2_signature':dict(Counter(map(len,fo))),'all_F2_cusps_lift_to_same_global_degree':True}
(P/'data'/'cusp133.json').write_text(json.dumps(out,indent=2)+'\n');print('NONIC133_Z_CUSPS_PASS',out)
