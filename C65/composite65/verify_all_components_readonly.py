#!/usr/bin/env python3
"""Independent complete arithmetic for ordinary C65 degrees7/8/9 only."""
from pathlib import Path
import json,subprocess,sys,time,os,collections
H=Path(__file__).resolve().parent;t=time.monotonic();r=json.loads((H.parent/'data'/'composite65'/'canonical_R1_sieve_verified.json').read_text());assert r['status']=='PASS'and r['cases']==336 and not r['unresolved']and r['source_diamond_orbits']==7 and r['source_closed_points']==168
assert collections.Counter((a['positive_degree'],a['columns'],a['rank'])for a in r['results'])=={(49,288,121):96,(56,295,121):96,(63,302,121):144};assert{(a['source'],a['target_diamond_translate'],a['convention'])for a in r['results']}=={(s,k,c)for s in range(1,8)for k in range(24)for c in[1,2]};assert max(m for a in r['results']for fid,m in a['positive_support']+a['negative_support'])==3
commands=[['verify_inputs_readonly.py'],['verify_stage1_readonly.py'],['verify_stage2_readonly.py','--readonly'],['verify_hecke_transport_readonly.py','--readonly'],['verify_census_readonly.py','--readonly'],['verify_cusp_pole_degrees.py','--readonly'],['verify_cusp_action_readonly.py','--readonly'],['build_R1_product_tensor.py','--readonly'],['canonical_R1_sieve.py','--readonly'],['verify_R1_ranks_m4ri.py']]
env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1';env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
for args in commands:
 z=subprocess.run([sys.executable,str(H/args[0])]+args[1:],text=True,env=env,timeout=300);assert z.returncode==0,(args,z.stdout,z.stderr);print('CHECK_PASS',args[0],flush=True)
print('NONIC65_ORDINARY_COMPONENTS_READONLY_PASS cases336 rank121 scope7_8_9_only elapsed',time.monotonic()-t)
