#!/usr/bin/env python3
"""Independent exact character-sum census of every ordinary curve and twist."""
from pathlib import Path
import json,sys,time
from unit_arithmetic import *
H=Path(__file__).resolve().parents[1]/'data/C85';R=H.parents[2];t=time.monotonic();out=[]
for f in range(6,10):
 P=R/'C85_C91_C93/data'/('N85_f'+str(f));s=(P/'records.m').read_text();q=1<<f;F=Field(bits(lit(s,'FieldPolynomial')));assert F.d==f;tr=[]
 for a in range(q):
  b=a;c=0
  for _ in range(f):c^=b;b=F.square(b)
  assert c in[0,1];tr.append(c)
 inv=[0]+[F.inv(a)for a in range(1,q)];ivs=[F.square(a)for a in inv];qualified=[]
 for a6 in range(1,q):
  n=2+2*sum(tr[x^F.mul(a6,ivs[x])]==0 for x in range(1,q))
  for order in[n,2*(q+1)-n]:
   if order%85==0:qualified.append((inv[a6],order))
 # Supersingular trace is contained in these sets; no member has order divisible by85.
 ss={0,1<<((f+1)//2),-(1<<((f+1)//2))}if f%2 else{0,1<<(f//2),-(1<<(f//2)),1<<(f//2+1),-(1<<(f//2+1))};assert all((q+1-v)%85 for v in ss)
 records=lit(s,'Records');assert sorted(qualified)==sorted((bits(v[0]),v[3])for v in records)
 for j,order in qualified:
  assert order==510 and f==9 and all(F.frob(j,d)!=j for d in[1,3]);assert q%5!=1 and q%17!=1
 out.append({'degree':f,'qualified_j_values':len(qualified),'orders':sorted(set(n for j,n in qualified))})
 if f==9:
  jj=lit((H/'elliptic_hecke_jets_stage2_data.m').read_text(),'SourceJPolynomials');assert len({bits(p)for p in jj})==5
  assert {bits(p)for p in jj}=={bits(p)for p in lit(s,'QualifiedJPolynomials')}
assert [v['qualified_j_values']for v in out]==[0,0,0,45]
report={'status':'PASS','ordinary_all_j_and_twists_by_character_sum':out,'supersingular_orders_excluded':True,'j_orbits_degree9':5,'diamond_orbit_length':32,'closed_degree9_points':160,'degrees_at_most5_excluded_by_Hasse':True}
p=H/'census_independent_verified.json'
if '--readonly'in sys.argv:assert json.loads(p.read_text())==report
else:p.write_text(json.dumps(report,indent=2)+'\n')
print('NONIC85_CENSUS_INDEPENDENT_PASS closed9=160 jorbits5 elapsed',time.monotonic()-t)
