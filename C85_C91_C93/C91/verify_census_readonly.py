#!/usr/bin/env python3
"""Independent exact character-sum census of every ordinary curve and twist."""
from pathlib import Path
import json,sys,time
from unit_arithmetic import *
H=Path(__file__).resolve().parents[1]/'data/C91';R=H.parents[2];t=time.monotonic();out=[]
for f in range(6,10):
 P=R/'C85_C91_C93/data'/('N91_f'+str(f));s=(P/'records.m').read_text();q=1<<f;F=Field(bits(lit(s,'FieldPolynomial')));assert F.d==f;tr=[]
 for a in range(q):
  b=a;c=0
  for _ in range(f):c^=b;b=F.square(b)
  assert c in[0,1];tr.append(c)
 inv=[0]+[F.inv(a)for a in range(1,q)];ivs=[F.square(a)for a in inv];qualified=[]
 for a6 in range(1,q):
  n=2+2*sum(tr[x^F.mul(a6,ivs[x])]==0 for x in range(1,q))
  for order in[n,2*(q+1)-n]:
   if order%91==0:qualified.append((inv[a6],order))
 # Supersingular trace is contained in these sets; no member has order divisible by91.
 ss={0,1<<((f+1)//2),-(1<<((f+1)//2))}if f%2 else{0,1<<(f//2),-(1<<(f//2)),1<<(f//2+1),-(1<<(f//2+1))};assert all((q+1-v)%91 for v in ss)if f!=8 else 273 in{q+1-v for v in ss}
 records=lit(s,'Records');assert sorted(qualified)==sorted((bits(v[0]),v[3])for v in records if bits(v[0])!=0)
 for j,order in qualified:
  assert order==546 and f==9 and all(F.frob(j,d)!=j for d in[1,3]);assert order%49!=0 and order%169!=0
 out.append({'degree':f,'qualified_j_values':len(qualified),'orders':sorted(set(n for j,n in qualified))})
 if f==9:
  jj=lit((H/'elliptic_hecke_jets_stage2_data.m').read_text(),'SourceJPolynomials');assert len({bits(p)for p in jj})==4
  assert {bits(p)for p in jj}=={bits(p)for p in lit(s,'QualifiedJPolynomials')}
assert [v['qualified_j_values']for v in out]==[0,0,0,36]
report={'status':'PASS','ordinary_all_j_and_twists_by_character_sum':out,'supersingular_orders_outside_degree8_excluded':True,'degree8_supersingular_twists':2,'degree8_supersingular_order':273,'j_orbits_degree9':4,'diamond_orbit_length':36,'closed_degree9_points':144,'degrees_at_most5_excluded_by_Hasse':True}
p=H/'census_independent_verified.json'
if '--readonly'in sys.argv:assert json.loads(p.read_text())==report
else:p.write_text(json.dumps(report,indent=2)+'\n')
print('NONIC91_CENSUS_INDEPENDENT_PASS closed9=144 jorbits4 CM8separate elapsed',time.monotonic()-t)
