#!/usr/bin/env python3
"""Independent primitive-column cusp and Frobenius-orbit audit at level44."""
from pathlib import Path
from math import gcd
from collections import Counter
import json,sys
N=44;units=[u for u in range(N) if gcd(u,N)==1]
def key(a,c):
 d=gcd(c,N);a%=d
 return min((a,c%N),((-a)%d,(-c)%N))
cusps={key(a,c) for c in range(N) for a in range(gcd(c,N)) if gcd(a,gcd(c,N))==1}
assert len(cusps)==50
# Galois acts on the numerator. The cusp width is N/gcd(c,N).
def gal(P,u):return key(P[0]*u,P[1])
def orbits(U):
 todo=set(cusps);ans=[]
 while todo:
  P=min(todo);orb=frozenset(gal(P,u) for u in U);assert orb<=todo;todo-=orb;ans.append(orb)
 return ans
powers=[];q=1
while q not in powers:powers.append(q);q=q*3%N
G=orbits(units);F=orbits(powers);gs=Counter(map(len,G));fs=Counter(map(len,F))
assert gs=={1:15,2:5,5:1,10:2} and fs=={1:15,2:5,5:3,10:1}
small=[o for o in F if len(o)<=4];assert len(small)==20
assert all(o in G for o in small)
assert all(frozenset(gal(P,u) for P in o)==o for o in small for u in units)
rows=[]
for w in sorted({N//gcd(P[1],N) for P in cusps}):
 gg=[o for o in G if N//gcd(next(iter(o))[1],N)==w];ff=[o for o in F if N//gcd(next(iter(o))[1],N)==w]
 rows.append({'width':w,'global_degrees':sorted(map(len,gg)),'finite_degrees':sorted(map(len,ff))})
r={'level':44,'geometric_cusps':50,'global_signature':dict(sorted(gs.items())),'F3_signature':dict(sorted(fs.items())),'width_rows':rows,'all_F3_cusps_degree_le4_lift_as_rational_closed_divisors':True,'all_cyclotomic_operators_fix_each_relevant_finite_orbit':True,'complete':True}
if '--write' in sys.argv:(Path(__file__).resolve().parent/'data'/'cusp_audit44.json').write_text(json.dumps(r,indent=2)+'\n')
print('C44_CUSP_AUDIT_PASS '+json.dumps(r))
