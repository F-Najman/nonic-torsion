#!/usr/bin/env python3
"""Primitive-column global/Frobenius cusp audit for X1(35)."""
from pathlib import Path
from collections import Counter
from math import gcd
import json
P=Path(__file__).resolve().parents[2]/'C35/data/presentation_v2';n=35
labels={(a,b) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1}
def orbit(x):
 a,b=x;return frozenset((s*(a+k*b)%n,s*b%n) for s in(-1,1) for k in range(n))
classes=sorted(set(map(orbit,labels)),key=lambda c:sorted(c));assert len(classes)==48
index={x:i for i,c in enumerate(classes) for x in c};reps=[min(c) for c in classes];units=[u for u in range(n) if gcd(u,n)==1]
gal={u:[index[(u*a%n,b)] for a,b in reps] for u in units}
remaining=set(range(48));glob=[]
while remaining:
 i=min(remaining);o={gal[u][i] for u in units};remaining-=o;glob.append(frozenset(o))
remaining=set(range(48));frob=[]
while remaining:
 i=min(remaining);o=[]
 while i not in o:o.append(i);remaining.remove(i);i=gal[3][i]
 assert i==o[0];frob.append(frozenset(o))
assert set(glob)==set(frob)
assert Counter(map(len,glob))=={1:12,4:3,6:2,12:1}
D={'level':35,'geometric_cusps':48,'global_signature':dict(Counter(map(len,glob))),'F3_signature':dict(Counter(map(len,frob))),'global_and_F3_orbit_partitions_equal':True,'all_F3_cusps_lift_to_same_global_degree':True,'FB_cusp_degrees':[1]*12+[4]*3+[6]*2,'missing_FB_cusp_degree':12,'complete':True}
(P/'cusp35.json').write_text(json.dumps(D,indent=2)+'\n');print('NONIC35_CUSP_AUDIT_PASS '+json.dumps(D))
