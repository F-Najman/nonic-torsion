#!/usr/bin/env python3
"""Read-only complete C55 degree-nine certificate."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
H=Path(__file__).resolve().parent;B=H.parent
def native(tag,marker):
 p=H/'logs'/(tag+'.log');s=p.read_text()
 assert marker in s and not re.search('Runtime error|User error|Assertion failed|syntax error',s,re.I)
 return s
s=native('N55_good9_retry','NONIC55_GOOD9_DONE all_places=40 representatives=2 complete=true')
assert re.findall(r'NONIC55_GOOD9_CASE index=(\d+) principal=(\w+)',s)==[('1','false'),('2','false')]
assert 'NONIC55_FULL_NATIVE_EQUAL degree=36'in s
native('independent_counts55_retry','closed_place_counts=[ 0, 0, 0, 0, 0, 0, 20, 0, 40 ]')
s=native('N55_coupled_FI_RR_v2','NONIC55_COUPLED_FI_AND_EXCEPTION_RR_DONE tested280=true exceptions8_h0one=true')
assert 'RANK0_CARTIER_DONE N=55 D=9 tested=280 total=280 exceptions=8 complete=true'in s
assert len(re.findall('NONIC55_COUPLED_EXCEPTION_RR vector=.*?h0=1',s))==8
native('N55_quotient_injectivity_D9_v2','NONIC55_QUOTIENT_INJECTIVITY_AND_RR_DONE')
G=B/'common';data=H/'data'/'N55_D9_coupled_v2_kernels.m'
outputs={}
for script in ['replay_cartier_projection.py','replay_rr_kernels.py','replay_fi_jets.py']:
 r=subprocess.run([sys.executable,str(G/script),str(data)],capture_output=True,text=True,check=True);outputs[script]=r.stdout;print(r.stdout.strip())
assert "'projected_jet_matrices': 38"in outputs['replay_cartier_projection.py']
r=json.loads(outputs['replay_fi_jets.py'].split('FI_REPLAY_DONE ')[1]);assert r['complete']and r['degree']==9 and r['representatives']==r['covered']==r['enumerated']==280 and r['exceptions']==8
for ln in outputs['replay_fi_jets.py'].splitlines():
 if ln.startswith('FI_REPLAY_EXCEPTION'):
  e=json.loads(ln.split('FI_REPLAY_EXCEPTION ')[1]);assert e['vector'][0]==0 and e['rank']==8
for script,marker in [('N55_tate_cusps.py','NONIC_ORDER5_QUOTIENT_CUSP_AUDIT_DONE 55'),('N55_class_replay.py','NONIC55_CUSPIDAL_IMAGE_HNF_REPLAY_PASS')]:
 r=subprocess.run([sys.executable,str(H/script)],capture_output=True,text=True,check=True);assert marker in r.stdout;print(marker)
s=(H/'logs'/'degree8'/'large_levels'/'C55_inputs'/'N55_D8_NC7_q3.log').read_text();assert 'principal=false'in s.replace('principal= false','principal=false')or'principal false'in s
print('NONIC55_READONLY_CERTIFICATE_PASS genus81 quotient17 good9_orbits2 all280_FI8exceptions_resolved rational2image32 complete=true')
