#!/usr/bin/env python3
"""Complete read-only C85 certificate; no Magma launch and no evidence mutation."""
from pathlib import Path
import hashlib,json,subprocess,sys,time,re
C=Path(__file__).resolve().parent;R=C.parent;H=C/'data/C85';K=C/'C85';t=time.monotonic();m=json.loads((H/'manifest.json').read_text());assert m['degree']==9 and m['group']==[1,85]and m['status']=='COMPLETE EXCLUSION'
supp=json.loads((H/'portability_supplement/manifest.json').read_text());assert supp['base_manifest_sha256']==hashlib.sha256((H/'manifest.json').read_bytes()).hexdigest()
for f in range(6,10):
 p=C/'data'/('N85_f'+str(f));assert'Complete:=true;'in(p/'records.m').read_text()
s=(C/'data/character_probes/N85_winding_data.m').read_text();assert'ZeroLFactors:=2;'in s and'FiniteProjectiveDiamondOrder:=8;'in s
r=json.loads((H/'canonical_R1_sieve_verified.json').read_text());assert r['status']=='PASS'and r['cases']==320 and not r['unresolved']and r['source_closed_points']==160;assert{(v['positive_degree'],v['columns'],v['rank'])for v in r['results']}=={(72,519,225)};assert max(c for v in r['results']for fid,c in v['positive_support']+v['negative_support'])==3
assert{(v['source'],v['target_diamond_translate'],v['convention'])for v in r['results']}=={(s,k,c)for s in range(1,6)for k in range(32)for c in[1,2]}
fi=json.loads((H/'FI_positive_independent_verified.json').read_text());assert fi['rank']==122 and fi['relation_kernel_dimension']==11 and fi['nonzero_kernel_vectors_exhausted']==2047 and fi['minimum_cusp_prefix_cost']==32 and fi['all_degree9_prefix_configurations_independent']
commands=[['verify_stage1_readonly.py'],['verify_stage2_readonly.py'],['verify_hecke_transport_readonly.py','--readonly'],['verify_census_readonly.py','--readonly'],['verify_cusp_pole_degrees.py','--readonly'],['verify_cusp_action_readonly.py','--readonly'],['portability_supplement/verify_FI_positive_portable_readonly.py','--readonly'],['build_R1_product_tensor.py','--readonly'],['canonical_R1_sieve.py','--readonly'],['verify_R1_ranks_m4ri.py']]
for args in commands:
 z=subprocess.run([sys.executable,str(K/args[0])]+args[1:],capture_output=True,text=True,timeout=300);assert z.returncode==0,(args,z.stdout,z.stderr);print('CHECK_PASS',args[0],flush=True)
print('NONIC85_COMPLETE_PORTABLE_READONLY_CERTIFICATE_PASS',json.dumps({'files':len(m['files'])+len(supp['files']),'selected_native_exits_zero':8,'binary_matrices':320,'matrix_shape':[225,519],'matrix_rank':225,'allbad_normalized_configurations':19070333,'genuine_ordinary_torsion_families':20,'claim':'No point of order85 over any degree9 number field','manifest_sha256':hashlib.sha256((H/'manifest.json').read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t},sort_keys=True))
