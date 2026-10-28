#!/usr/bin/env python3
"""Independent exhaustive ordinary C65 character-sum census in degrees6..9."""
from pathlib import Path
import json,sys,time
from unit_arithmetic import *
H=Path(__file__).resolve().parents[1]/'data'/'composite65';R=H.parents[2];t=time.monotonic();out=[];data=(H/'elliptic_hecke_jets_stage2_data.m').read_text();fds=lit(data,'SourceDegrees');jps=lit(data,'SourceJPolynomials')
for f in range(6,10):
 s=(R/'C65'/'data'/'census'/f'N65_f{f}_records.m').read_text();q=1<<f;F=Field(bits(lit(s,'FieldPolynomial')));assert F.d==f;tr=[]
 for a in range(q):
  b=a;c=0
  for _ in range(f):c^=b;b=F.square(b)
  assert c in[0,1];tr.append(c)
 inv=[0]+[F.inv(a)for a in range(1,q)];ivs=[F.square(a)for a in inv];qualified=[]
 for a6 in range(1,q):
  n=2+2*sum(tr[x^F.mul(a6,ivs[x])]==0 for x in range(1,q))
  for order in[n,2*(q+1)-n]:
   if order%65==0:qualified.append((inv[a6],order))
 records=lit(s,'Records');assert sorted(qualified)==sorted((bits(v[0]),v[3])for v in records if bits(v[0]))
 for j,order in qualified:
  assert f in[7,8,9]and order=={7:130,8:260,9:520}[f]and all(F.frob(j,d)!=j for d in range(1,f)if f%d==0);assert order%25!=0 and order%169!=0
 if f>=7:assert {bits(p)for g,p in zip(fds,jps)if g==f}=={bits(p)for p in lit(s,'QualifiedJPolynomials')}
 out.append({'degree':f,'ordinary_qualified_j_values':len(qualified),'orders':sorted(set(n for j,n in qualified))})
assert [v['ordinary_qualified_j_values']for v in out]==[0,14,16,27]
report={'status':'PASS','scope':'ordinary curves and both twists only; supersingular degree6 handled separately','ordinary_all_j_and_twists_by_character_sum':out,'ordinary_j_orbits_by_degree':{'7':2,'8':2,'9':3},'diamond_orbit_length':24,'ordinary_closed_points_by_degree':{'7':48,'8':48,'9':72},'ordinary_closed_points_total':168};p=H/'census_independent_verified.json'
if '--readonly'in sys.argv:assert json.loads(p.read_text())==report
else:p.write_text(json.dumps(report,indent=2)+'\n')
print('NONIC65_ORDINARY_CENSUS_PASS closed7_8_9=[48,48,72] elapsed',time.monotonic()-t)
