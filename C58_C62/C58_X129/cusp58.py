#!/usr/bin/env python3
"""Primitive-column global/Frobenius cusp audit for X1(29)."""
from pathlib import Path
from collections import Counter
from math import gcd
import json
P=Path(__file__).resolve().parent;n=29
labels={(a,b) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1}
def orbit(x):
 a,b=x;return frozenset((s*(a+k*b)%n,s*b%n) for s in(-1,1) for k in range(n))
classes=sorted(set(map(orbit,labels)),key=lambda c:sorted(c));assert len(classes)==28
index={x:i for i,c in enumerate(classes) for x in c};reps=[min(c) for c in classes];units=[u for u in range(n) if gcd(u,n)==1]
gal={u:[index[(u*a%n,b)] for a,b in reps] for u in units}
remaining=set(range(28));glob=[]
while remaining:
 i=min(remaining);o={gal[u][i] for u in units};remaining-=o;glob.append(frozenset(o))
remaining=set(range(28));frob=[]
while remaining:
 i=min(remaining);o=[]
 while i not in o:o.append(i);remaining.remove(i);i=gal[3][i]
 assert i==o[0];frob.append(frozenset(o))
assert set(glob)==set(frob)
assert Counter(map(len,glob))=={1:14,14:1}
D={'level':29,'geometric_cusps':28,'global_signature':dict(Counter(map(len,glob))),'F3_signature':dict(Counter(map(len,frob))),'global_and_F3_orbit_partitions_equal':True,'all_F3_cusps_lift_to_same_global_degree':True,'FB_cusp_degrees':[1]*14,'missing_FB_cusp_degree':14,'complete':True}
(P.parent/'data'/'C58_X129'/'cusp58.json').write_text(json.dumps(D,indent=2)+'\n');print('NONIC58_CUSP_AUDIT_PASS '+json.dumps(D))
