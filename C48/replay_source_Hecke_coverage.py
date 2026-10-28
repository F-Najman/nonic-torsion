#!/usr/bin/env python3
"""Read-only exact source-case and cusp-complement binding for degree9."""
from pathlib import Path
import ast,json,runpy,re
from math import gcd
from collections import Counter
P=Path(__file__).resolve().parent

def literal(path):return {n.targets[0].id:ast.literal_eval(n.value) for n in ast.parse(path.read_text().replace('[*','[').replace('*]',']')).body}
def magma(path):
 out={}
 for stmt in path.read_text().split(';'):
  if ':=' not in stmt:continue
  name,value=stmt.split(':=',1);value=value.strip().replace('[*','[').replace('*]',']');value={'true':'True','false':'False'}.get(value,value);out[name.strip()]=ast.literal_eval(value)
 return out
D=json.loads((P/'data/first_D9_filter.json').read_text());S=runpy.run_path(str(P/'data/degree8/source_partition_inputs_complete.py'));L=runpy.run_path(str(P/'data/degree8/clean_run/lower_support_arithmetic.py'))
expected13=[]
for ds,ix in D['low_survivors']:expected13.append([ds,ix,[S['representative_pairs'][ds[0]][ix[0]-1],L['all_source_pairs'][ds[1]][ix[1]-1]]])
for d,i,cid in D['cached_high_survivors']:
 if d==8:expected13.append([[8],[i],[S['representative_pairs'][8][i-1]]])
assert len(expected13)==26 and json.loads((P/'data/all_hecke13_cases.json').read_text())['cases']==expected13
for folder,ids in [('hecke13_0_v2',list(range(1,14))),('hecke13_1',list(range(14,27)))]:
 I=magma(P/'data'/folder/'inputs.m');O=literal(P/'data'/folder/'results.py');assert I['CaseIDs']==ids and I['Cases']==[expected13[i-1] for i in ids]
 assert I['FieldPolynomials']==[S['field_polynomials'][d] for d in range(3,9)]
 assert O['q']==13 and O['complete'] and O['remaining']==[] and len(O['tests'])==13
 for row,i in zip(O['tests'],ids):assert row==[i,expected13[i-1][0],expected13[i-1][1],0]
P9=magma(P/'data/source9_survivors/records.m');assert P9['Complete'] and len(P9['Source9Cases'])==24
passindices=json.loads((P/'data/class9_survivors.json').read_text())['target_Aut16_orbit_indices'];assert magma(P/'data/source9_survivors/inputs.m')['PassingIndices']==passindices
assert len({tuple(r[1]+r[2]) for r in P9['Source9Cases']})==24 and sum(r[3] for r in P9['Source9Cases'])==192
for folder,ids in [('hecke7_0',list(range(1,13))),('hecke7_1',list(range(13,25)))]:
 I=magma(P/'data'/folder/'inputs.m');O=literal(P/'data'/folder/'results.py');assert I['SourceFieldPolynomial']==P9['Source9FieldPolynomial'];assert I['CaseIDs']==ids and I['Cases']==[P9['Source9Cases'][i-1] for i in ids]
 assert O['q']==7 and O['complete'] and O['remaining']==[] and len(O['tests'])==12
 for row,i in zip(O['tests'],ids):assert row==[i,P9['Source9Cases'][i-1][0],0]
old4=[47,85,171,243]
for d in [6,7]:
 wanted=[i for dd,i,cid in D['cached_high_survivors'] if dd==d and (d==7 or i not in old4)];assert len(wanted)==(10 if d==6 else 8)
 I=magma(P/f'data/hecke29_d{d}/inputs.m');O=literal(P/f'data/hecke29_d{d}/results.py');ids=list(range(1,len(wanted)+1))
 assert I['SourceFieldPolynomial']==S['field_polynomials'][d] and I['CaseIDs']==ids and I['Cases']==[[i,*S['representative_pairs'][d][i-1]] for i in wanted]
 assert O['q']==29 and O['complete'] and O['remaining']==[] and len(O['tests'])==len(wanted)
 for row,i,sourceidx in zip(O['tests'],ids,wanted):assert row==[i,sourceidx,0]
assert set(old4).issubset({i for d,i,cid in D['cached_high_survivors'] if d==6})
for sh,idx in enumerate(old4):
 s=(P/f'logs/degree8/high_degree6_shard{sh}.log').read_text();assert f'C48_D6_HEC_DONE source_index={idx} q=29 principal=false complete=true' in s
# Independent Galois/Frobenius audit of every finite cusp orbit.
pts=[]
for c in range(25):
 d=gcd(48,c)
 for a in range(1,d+1):
  if gcd(a,d)!=1 or (c in(0,24) and a>d-a):continue
  pts.append((c,a%d))
def gal(pt,q):
 c,a=pt;d=gcd(48,c);a=q*a%d
 if c in (0,24):a=min(a,(-a)%d)
 return c,a
seen=set();orb=[]
for pt in pts:
 if pt in seen:continue
 o=[];t=pt
 while t not in o:o.append(t);seen.add(t);t=gal(t,5)
 orb.append(set(o))
assert len(pts)==56 and Counter(map(len,orb))=={1:16,2:12,4:4}
C=json.loads((P/'data/fresh_cusp_Hecke_audit.json').read_text());assert C['geometric_cusps']==[list(v) for v in pts] and C['Frobenius5_orbits']==[[list(v) for v in sorted(o)] for o in orb]
for q in (7,13,29):
 bad=[sorted(o) for o in orb if {gal(v,q) for v in o}!=o]
 assert C['operators'][str(q)]['nonstable_cusp_orbits']==[[list(v) for v in o] for o in bad]
 if q==29:assert not bad and pow(5,3,48)==29
 if q==13:assert Counter(map(len,bad))=={2:2,4:2}
 if q==7:assert Counter(map(len,bad))=={1:4,2:6,4:4}
assert 26+24+10+8+4==72
print('NONIC48_SOURCE_COVERAGE_REPLAY_PASS '+json.dumps({'fresh_direct_RR_tests':68,'retained_valid_A29_tests':4,'all72_source_cases_excluded':True,'q29_all_cusps_fixed':True,'q13_degree1_cusps_fixed':True,'q13_degree2_not_assumed_fixed':True,'complete':True}))
