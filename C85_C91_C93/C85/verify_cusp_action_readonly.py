#!/usr/bin/env python3
"""Exact squarefree-level cusp/Frobenius/involution arithmetic."""
from pathlib import Path
import math,json,sys
H=Path(__file__).resolve().parents[1]/'data/C85';N=85
# A cusp is (a,c), c modulo N, a modulo gcd(c,N), modulo simultaneous sign.
def normal(a,c):
 c%=N;g=math.gcd(c,N);a%=g;assert math.gcd(a,g)==1
 return min((a,c),((-a)%g,(-c)%N))
cusps={normal(a,c)for c in range(N)for a in range(math.gcd(c,N))if math.gcd(a,math.gcd(c,N))==1};assert len(cusps)==128
def frob(v):a,c=v;return normal(2*a,c)
def h(v):a,c=v;return normal(16*a,16*c)
assert all(h(h(v))==v for v in cusps)
seen=set();orbits=[]
for v in sorted(cusps):
 if v in seen:continue
 o=[];w=v
 while w not in o:o.append(w);seen.add(w);w=frob(w)
 assert w==v;orbits.append(o)
assert sorted(map(len,orbits))==[1]*32+[4]*8+[8]*8
rat={o[0]for o in orbits if len(o)==1};assert all(h(v)!=v and h(v)in rat for v in rat)
assert all({h(v)for v in o}==set(o)for o in orbits if len(o)>1)
assert all((N//math.gcd(v[1],N)==85)==(v in rat)for v in cusps)
out={'status':'PASS','level':85,'geometric_cusps':128,'Frobenius_closed_degrees':{'1':32,'4':8,'8':8},'h_scalar':16,'h_order':2,'h_free_on_rational_cusps':True,'h_fixes_each_nonrational_closed_cusp':True,'all_cusp_pairs':[list(v)for v in sorted(cusps)],'Frobenius_orbits':[[list(v)for v in o]for o in orbits]}
p=H/'cusp_action_verified.json'
if '--readonly'in sys.argv:assert json.loads(p.read_text())==out
else:p.write_text(json.dumps(out,indent=2)+'\n')
print('NONIC85_CUSP_ACTION_PASS degrees1^32_4^8_8^8 h16_free_rat_fixes_nonrat=true')
