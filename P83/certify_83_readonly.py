#!/usr/bin/env python3
"""Read-only complete C83 packet replay; no native launch."""
from pathlib import Path
import hashlib,json,subprocess,sys,time,re,ast,gzip
C=Path(__file__).resolve().parent;ROOT=C.parent;P=C/'data';Q=C/'canonical_units';QD=P/'canonical_units';QL=C/'logs/canonical_units';G=ROOT/'common/logs';t=time.monotonic();m=json.loads((P/'manifest.json').read_text());assert m['degree']==9 and m['group']==[1,83]and m['status']=='COMPLETE EXCLUSION'
assert'NONIC83_CANONICAL_UNIT_STAGE1_COMPLETE generic_rank246=true no_function_field_normalization=true'in(QL/'canonical_unit_stage1.log').read_text()
log=(QL/'elliptic_hecke_jets_stage2_v2.log').read_text();assert'NONIC83_ELLIPTIC_HECKE_JETS_STAGE2_COMPLETE'in log;assert len(re.findall(r'target_closed_degrees=\[ 9, 9, 18 \] norm_multiplicities=\[ 1, 1, 1 \]',log))==5
log=(G/'finite_remaining_census.log').read_text()
for d in[7,8,9]:
 line=next(x for x in log.splitlines()if f'N=83 degree={d} 'in x);assert'complete=true'in line
 if d in[7,8]:assert'closed_places=0'in line
 else:assert'closed_places=205'in line and ast.literal_eval(re.search(r'diamond_orbit_sizes=(\[.*?\])',line)[1])==[41]*5 and ast.literal_eval(re.search(r'j_degrees=(\[.*?\])',line)[1])==[9]*5
assert re.search(r'TORSION n=83 p=2.*= \[ 9 \]',(G/'finite_remaining_primes_49.log').read_text())
fi=json.loads((P/'FI_verified.json').read_text());assert fi['p']==83 and fi['H_index']==41 and fi['full_binary_rank']==169 and fi['relation_kernel_dimension']==0 and fi['prefix_configurations']==134858144
r=json.loads((QD/'canonical_R1_sieve_verified.json').read_text());assert r['status']=='PASS'and r['cases']==410 and not r['unresolved']and r['source_closed_points']==205
assert {(v['positive_degree'],v['columns'],v['rank'])for v in r['results']}=={(72,563,247)};assert max(mult for v in r['results']for fid,mult in v['positive_support']+v['negative_support'])==3
assert {(v['source'],v['target_diamond_translate'],v['convention'])for v in r['results']}=={(s,k,c)for s in range(1,6)for k in range(41)for c in[1,2]}
commands=[[sys.executable,str(Q/p)]for p in['verify_stage1_readonly.py','verify_stage2_readonly.py','verify_R1_ranks_m4ri.py']]
commands += [[sys.executable,str(Q/'build_R1_product_tensor.py'),'--readonly'],[sys.executable,str(Q/'canonical_R1_sieve.py'),'--readonly'],[sys.executable,str(ROOT/'common/verify_FI_relations.py'),str(P/'FI_rows.txt')]]
for cmd in commands:
 z=subprocess.run(cmd,capture_output=True,text=True,timeout=300);assert z.returncode==0,(cmd,z.stdout,z.stderr);print('CHECK_PASS',Path(cmd[1]).name,flush=True)
print('NONIC83_COMPLETE_READONLY_CERTIFICATE_PASS',json.dumps({'files':len(m['files']),'selected_native_exits_zero':5,'binary_matrices':410,'matrix_shape':[247,563],'matrix_rank':247,'allbad_configurations':134858144,'genuine_ordinary_torsion_families':20,'claim':'No point of order83 over any degree9 number field','manifest_sha256':hashlib.sha256((P/'manifest.json').read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t},sort_keys=True))
