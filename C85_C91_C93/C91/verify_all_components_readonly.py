#!/usr/bin/env python3
"""Complete independent C91 arithmetic replay, with no Magma launch."""
from pathlib import Path
import sys,json,subprocess,time,os
C=Path(__file__).resolve().parent;H=C.parent/'data/C91';t=time.monotonic();r=json.loads((H/'canonical_R1_sieve_verified.json').read_text());assert r['status']=='PASS'and r['cases']==288 and not r['unresolved']and r['source_diamond_orbits']==4 and r['source_closed_points']==144
assert{(v['positive_degree'],v['columns'],v['rank'])for v in r['results']}=={(63,590,265)};assert{(v['source'],v['target_diamond_translate'],v['convention'])for v in r['results']}=={(s,k,c)for s in range(1,5)for k in range(36)for c in[1,2]}
assert max(c for v in r['results']for fid,c in v['positive_support']+v['negative_support'])==3
c=json.loads((H/'CM8_action_verified.json').read_text());assert c['h_fixes_all_CM8_closed_points']and c['closed_points']==3 and c['Tate_normal_form_inverse_regular_and_exact_on_all_retained_points']
commands=[['verify_inputs_readonly.py'],['verify_stage1_readonly.py'],['verify_stage2_readonly.py'],['verify_hecke_transport_readonly.py','--readonly'],['verify_census_readonly.py','--readonly'],['verify_CM8_action.py','--readonly'],['verify_cusp_pole_degrees.py','--readonly'],['verify_cusp_action_readonly.py','--readonly'],['build_R1_product_tensor.py','--readonly'],['canonical_R1_sieve.py','--readonly'],['verify_R1_ranks_m4ri.py']]
env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1';env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
for args in commands:
 z=subprocess.run([sys.executable,str(C/args[0])]+args[1:],capture_output=True,text=True,timeout=300,env=env);assert z.returncode==0,(args,z.stdout,z.stderr);print('CHECK_PASS',args[0],flush=True)
print('NONIC91_ALL_COMPONENTS_READONLY_PASS cases288 rank265 CM8h27=true elapsed',time.monotonic()-t)
