#!/usr/bin/env python3
"""Read-only complete C71 packet replay; no Magma job is launched."""
from pathlib import Path
import hashlib,json,subprocess,sys,time,re,ast
ROOT=Path(__file__).resolve().parents[1];C=ROOT/'common';P=ROOT/'P71';Q=P/'quotient26';D=P/'data';QD=D/'quotient26';LG=P/'logs';t=time.monotonic()
log=(LG/'quotient26'/'native_control.log').read_text();assert'NONIC71_SIGMA_NATIVE_CONTROL order5=true fixes_j=true'in log and'NONIC71_RR_NATIVE_CONTROL_COMPLETE degree7 convention1 dimension0'in log and'NONIC71_RR_NATIVE_CONTROL_COMPLETE degree8 convention2 dimension0'in log and'NONIC71_NATIVE_CONTROLS_COMPLETE'in log
log=(C/'logs'/'finite_remaining_census.log').read_text()
for d in[7,8,9]:
 line=next(x for x in log.splitlines()if f'N=71 degree={d} 'in x);assert'complete=true'in line
 if d==9:assert'closed_places=0'in line
 else:
  orbits=ast.literal_eval(re.search(r'diamond_orbit_sizes=(\[.*?\])',line)[1]);js=ast.literal_eval(re.search(r'j_degrees=(\[.*?\])',line)[1]);assert sorted(zip(orbits,js))==([(5,1),(35,7)]if d==7 else[(35,8)])
log=(C/'logs'/'finite_remaining_primes_49.log').read_text();assert re.search(r'TORSION n=71 p=2.*= \[ 7, 8 \]',log)
fi=json.loads((D/'FI_relation_verified.json').read_text());assert fi['p']==71 and fi['H_index']==5 and fi['full_binary_rank']==25 and fi['relation_kernel_dimension']==0 and fi['prefix_configurations']==170
r=json.loads((QD/'canonical_ring_sieve_all_verified.json').read_text());assert r['status']=='PASS'and r['cases']==60 and not r['unresolved'];assert sorted({(a['positive_degree'],a['columns'],a['rank'])for a in r['results']})==[(28,128,75),(32,132,75)]
commands=[[sys.executable,str(Q/p)]for p in['verify_q_identities.py','verify_finite_keys.py','verify_canonical_ring_tensor_readonly.py','verify_ranks_m4ri.py']]
commands += [[sys.executable,str(Q/'canonical_ring_sieve.py'),'--readonly'],[sys.executable,str(C/'verify_FI_relations.py'),str(D/'FI_rows.txt')],[sys.executable,str(P/'certify_special7_readonly.py')]]
for cmd in commands:
 r=subprocess.run(cmd,capture_output=True,text=True,timeout=180);assert r.returncode==0,(cmd,r.stdout,r.stderr);print('CHECK_PASS',Path(cmd[1]).name,flush=True)
print('NONIC71_COMPLETE_READONLY_CERTIFICATE_PASS',json.dumps({'binary_matrices':60,'matrix_rank':75,'matrix_shapes':[[75,128],[75,132]],'allbad_configurations':170,'claim':'No point of order71 over any degree9 number field','elapsed_seconds':time.monotonic()-t},sort_keys=True))
