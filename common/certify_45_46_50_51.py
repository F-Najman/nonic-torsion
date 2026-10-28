#!/usr/bin/env python3
"""Read-only replay of four completed nonic cyclic exclusions."""
from pathlib import Path
import json,hashlib,re,subprocess,sys,ast
B=Path(__file__).resolve().parents[1]
LOG={'independent_counts_46_50_51':'common/logs/independent_counts_46_50_51.log','native_controls_short_fields':'common/logs/native_controls_short_fields.log','N51_good9_and_pure9':'C51/logs/N51_good9_and_pure9.log','N45_good9_and_pure9_retry':'C45/logs/N45_good9_and_pure9_retry.log','independent_counts45':'C45/logs/independent_counts45.log','N45_quotient_injectivity_control':'C45/logs/N45_quotient_injectivity_control.log'}
def native(tag,marker):
 p=B/LOG[tag];s=p.read_text()
 assert marker in s and not re.search(r'Runtime error|User error|Assertion failed|syntax error|failed memory request',s,re.I)
 return s
native('independent_counts_46_50_51','DEGREE9_FAST_COUNTS_DONE complete=true')
s=native('native_controls_short_fields','NONIC_NATIVE_CONTROLS_COMPLETE levels=46,50,51 source_degrees=9,6,9')
assert all(f'NONIC_NATIVE_HECKE_EQUAL N={n}'in s for n in[46,50,51])
s=native('N51_good9_and_pure9','NONIC51_GOOD9_AND_PURE9_DONE all_places=80 representatives=5 complete=true')
assert 'NONIC51_PURE9_CUSP_CONTROL rank=9 covered=16 dimensionU=62'in s
assert [(int(i),q,f)for i,q,f in re.findall(r'NONIC51_GOOD9_CASE index=(\d+) q=(\d+) principal=(\w+)',s)]==[(i,'5','false')for i in range(1,6)]
s=native('N45_good9_and_pure9_retry','NONIC45_GOOD9_AND_PURE9_DONE all_places=60 representatives=5 complete=true')
assert 'NONIC45_PURE9_CUSP_CONTROL rank=9 covered=12 dimensionU=39'in s and 'NONIC45_FULL_NATIVE_EQUAL degree=72'in s
assert [(int(i),q,f)for i,q,f in re.findall(r'NONIC45_GOOD9_CASE index=(\d+) q=(\d+) principal=(\w+)',s)]==[(i,'7','false')for i in range(1,6)]
native('independent_counts45','DEGREE9_FAST_COUNTS_DONE complete=true')
s=native('N45_quotient_injectivity_control','INJECTIVE_REDUCTION_3_15_DONE')
assert 'EQUAL_TWO_PRIMARY_ORDERS 8192 8192'in s
for name in ['N45','N51']:
 case='C45'if name=='N45'else'C51'
 data=B/(case+'/data/degree8/'+name+('_completion_jets.m'if name=='N45'else'_frozen_jets.m'))
 cmd=[sys.executable,str(B/'common/replay_cartier_projection.py'),str(data)]
 r=subprocess.run(cmd,capture_output=True,text=True,check=True);print(r.stdout.strip())
 r=subprocess.run([sys.executable,str(B/case/(name+'_replay_degree9.py')),'--readonly'],capture_output=True,text=True,check=True)
 assert name+'_DEGREE9_CUSP_REPLAY_DONE'in r.stdout;print(r.stdout.splitlines()[-1])
 j=json.loads((B/case/'data'/(name+'_degree9_FI_replay.json')).read_text())
 if name=='N51':assert j['tested_divisors']==3333248 and j['exceptions_count']==0
 else:
  assert j['tested_divisors']==377348 and j['exceptions_count']==8
  weights=ast.literal_eval(re.search(r'ds:=(.*?);',data.read_text())[1])
  assert all(e['rank']==8 and sum(a*w for a,w in zip(e['vector'],weights))==9 and all(w<=2 for a,w in zip(e['vector'],weights)if a)for e in j['exceptions'])
r=subprocess.run([sys.executable,str(B/'C45/N45_quotient_cusp_labels.py')],capture_output=True,text=True,check=True);assert 'NONIC45_QUOTIENT_RATIONAL_CUSPS_IDENTIFIED_PASS'in r.stdout
r=subprocess.run([sys.executable,str(B/'C46_C50/audit_46_50_D9.py')],capture_output=True,text=True,check=True);assert 'NONIC46_50_READONLY_COVERAGE_PASS'in r.stdout;print(r.stdout.strip())
s=(B/'C45/logs/degree8/N45_p2_good_only.log').read_text();assert 'RESULT type=[ 8 ] total=5 survivors=0'in s and 'total_survivors=0'in s and 'CC complement degree <= 1:'in s
print('NONIC45_46_50_51_READONLY_CERTIFICATE_PASS complete=true orders=45,46,50,51')
