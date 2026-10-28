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
index={x:i for i,c in enumerate(classes) for x in c};reps=[min(c) for c in classes]
def quotient(gen,order):
 hs={pow(gen,j,n) for j in range(order)};assert len(hs)==order
 qos=sorted({tuple(sorted({index[(h*a%n,h*b%n)] for h in hs})) for a,b in reps})
 assert all(len(o)==order for o in qos)
 qi={i:j for j,o in enumerate(qos) for i in o}
 fr=[qi[index[(2*reps[o[0]][0]%n,reps[o[0]][1])]] for o in qos]
 widths=[n//gcd(reps[o[0]][1],n) for o in qos]
 left=set(range(len(qos)));fos=[]
 while left:
  i=min(left);o=[]
  while i not in o:o.append(i);left.remove(i);i=fr[i]
  assert i==o[0] and len({widths[j] for j in o})==1
  fos.append(o)
 return qos,qi,fos,widths
Y,yqi,yfo,yw=quotient(64,3);Z,zqi,zfo,zw=quotient(78,18)
assert Counter(map(len,yfo))=={1:18,3:6,6:6}
zbyq={q:j for j,o in enumerate(zfo) for q in o}
admissible=set();image_info=[]
for o in yfo:
 if len(o)>3:continue
 image={zqi[Y[j][0]] for j in o};zj=zbyq[next(iter(image))]
 assert image==set(zfo[zj])
 admissible.add(zj);image_info.append([len(o),yw[o[0]],len(zfo[zj]),zw[zfo[zj][0]]])
types=Counter((len(zfo[i]),zw[zfo[i][0]]) for i in admissible)
assert types=={(1,133):3,(3,19):1}
alltypes=Counter((len(o),zw[o[0]]) for o in zfo)
assert alltypes=={(1,133):3,(1,7):3,(3,19):1,(3,1):1}
out={'level':133,'Y_diamond':64,'Z_diamond':78,'Y_F2_signature':dict(Counter(map(len,yfo))),'Z_F2_degree_width_counts':[[d,w,c] for (d,w),c in sorted(alltypes.items())],'admissible_Z_degree_width_counts':[[d,w,c] for (d,w),c in sorted(types.items())],'Y_low_degree_image_rows':image_info,'admissible_degree3_divisors':11,'complete':True}
(P/'data'/'cusp_image133.json').write_text(json.dumps(out,indent=2)+'\n');print('NONIC133_ADMISSIBLE_CUSP_IMAGE_PASS',out)
