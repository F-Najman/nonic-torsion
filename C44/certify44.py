#!/usr/bin/env python3
"""Read-only audit of the complete nonic C44 proof; no Magma launches."""
from pathlib import Path
import json,hashlib,subprocess,sys,re
P=Path(__file__).resolve().parent;H=P.parent
for j in [{'log':'C44/logs/N44_source_counts.log','completion_marker':'NONIC44_SOURCE_FROBENIUS_COUNTS_PASS'},{'log':'C44/logs/C44_canonical_build_v2.log','completion_marker':'NONIC44_CANONICAL_MODEL_CERTIFIED genus36=true degree70=true quadrics561=true'},{'log':'C44/logs/C44_finite_jets_all_v2.log','completion_marker':'NONIC44_FINITE_JETS_COMPLETE cases=226 source_distinct_keys=2260 all226=true'},{'log':'C44/logs/C44_native_control.log','completion_marker':'NONIC44_NATIVE_RR_CONTROL_DONE dimension=0 no_class_cache=true'},{'log':'C44/logs/C44_tensor22.log','completion_marker':'NONIC44_RING22_COMPLETE R2=105 R4=245 products11025=true'},{'log':'C44/logs/C44_matrix22_all.log','completion_marker':'NONIC44_CANONICAL22_COMPLETE cases=226 survivors=[] all_rank105=true'}]:
 text=(H/j['log']).read_text();assert j['completion_marker'] in text
 assert not re.search(r'Runtime error|User error|Assertion failed|Segmentation fault|Internal error',text)
for script,marker in [('audit_cusps44.py','C44_CUSP_AUDIT_PASS'),('audit_frobenius44.py','C44_FROBENIUS_AUDIT_PASS')]:
 r=subprocess.run([sys.executable,str(P/script)],capture_output=True,text=True,check=True);assert marker in r.stdout
O=json.loads((P/'data/degree8/N44_audit.json').read_text());assert O['degree']==8 and O['level']==44 and O['survivors']==0 and O['complete']
assert O['tested_divisors']==138 and [r['total'] for r in O['nonempty_types']]==[6,16,24,92]
def stream(path):
 text=path.read_text();dec=json.JSONDecoder();i=0;out=[]
 while i<len(text):
  while i<len(text) and text[i].isspace():i+=1
  if i==len(text):break
  x,i=dec.raw_decode(text,i);out.append(x)
 return out
J=stream(P/'data/sieve_all_inputs.normalized.jsonl');assert J[0]['source_places9']==2260 and J[0]['source_diamond_orbits']==226 and J[-1]['cases']==226 and J[-1]['source_distinct_keys']==2260 and J[-1]['all226']
assert sorted(r['case'] for r in J if r['type']=='source')==list(range(1,227))
for k in range(1,227):
 for side in ('T','U'):assert sum(r['degree']*r['multiplicity'] for r in J if r['type']=='jet' and r['case']==k and r['side']==side)==54
C=stream(P/'data/matrix22_all.raw.jsonseq');assert C[-1]['type']=='complete' and C[-1]['cases']==226 and C[-1]['survivors']==[]
R=[r for r in C if r['type']=='case'];assert [r['case'] for r in R]==list(range(1,227)) and all(r['rank']==105 for r in R)
for r in R:
 for k,n in [('F',105),('H',105),('ET2',105*54),('EU2',105*54),('quotient_projection',245*140),('matrix',105*194)]:assert len(r[k])==n and set(r[k])<={0,1,2}
A=json.loads((P/'data/matrix22_all_independent.json').read_text());assert A['complete'] and A['cases']==226 and A['all_obstruction_ranks']==105 and A['all_local_restrictions_reconstructed'] and A['all_product_annihilator_identities']
L=(P/'logs/matrix22_all_independent.log').read_text();assert L.count('C44_SAGE_MATRIX_CASE case=')==226 and 'C44_SAGE_MATRIX_REPLAY_PASS' in L
assert 'NONIC44_UNTOUCHED_NATIVE5_CONTROL_PASS all6=true' in (P/'logs/C44_native_control.log').read_text()
assert 'NONIC44_NATIVE_RR_CONTROL_DONE dimension=0 no_class_cache=true' in (P/'logs/C44_native_control.log').read_text()
if '--deep' in sys.argv:
 subprocess.run([sys.executable,str(P/'replay_matrices22_sage.py'),str(P/'data/matrix22_all.raw.jsonseq'),'--jets',str(P/'data/sieve_all_inputs.normalized.jsonl')],check=True)
print('NONIC_C44_READONLY_CERTIFICATE_PASS degree=9 inherited138_complement4=true source_places=2260 cases=226 all_rank105=true native_jobs=6 complete=true')
