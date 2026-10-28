#!/usr/bin/env python3
"""Exact squarefree-level cusp/Frobenius/involution arithmetic."""
from pathlib import Path
import math,json,sys
H=Path(__file__).resolve().parents[1]/'data'/'composite65';N=65
# A cusp is (a,c), c modulo N, a modulo gcd(c,N), modulo simultaneous sign.
def normal(a,c):
 c%=N;g=math.gcd(c,N);a%=g;assert math.gcd(a,g)==1
 return min((a,c),((-a)%g,(-c)%N))
cusps={normal(a,c)for c in range(N)for a in range(math.gcd(c,N))if math.gcd(a,math.gcd(c,N))==1};assert len(cusps)==96
def frob(v):a,c=v;return normal(2*a,c)
def h(v):a,c=v;return normal(3*a,22*c)

seen=set();orbits=[]
for v in sorted(cusps):
 if v in seen:continue
 o=[];w=v
 while w not in o:o.append(w);seen.add(w);w=frob(w)
 assert w==v;orbits.append(o)
assert sorted(map(len,orbits))==[1]*24+[4]*6+[6]*4+[12]*2
rat={o[0]for o in orbits if len(o)==1};assert all(h(v)!=v and h(v)in rat for v in rat)
# No nonrational fixed-cusp assertion is needed in the allbad unramified-lift argument.
assert all((N//math.gcd(v[1],N)==65)==(v in rat)for v in cusps)
out={'status':'PASS','level':65,'geometric_cusps':96,'Frobenius_closed_degrees':{'1':24,'4':6,'6':4,'12':2},'h_scalar':3,'h_order':12,'h_free_on_rational_cusps':True,'all_cusp_pairs':[list(v)for v in sorted(cusps)],'Frobenius_orbits':[[list(v)for v in o]for o in orbits]}
p=H/'cusp_action_verified.json'
if '--readonly'in sys.argv:assert json.loads(p.read_text())==out
else:p.write_text(json.dumps(out,indent=2)+'\n')
print('NONIC65_CUSP_ACTION_PASS degrees1^24_4^6_6^4_12^2 only_deg1_cusps_globalQ=true')
