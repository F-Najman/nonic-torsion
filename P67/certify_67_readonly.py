#!/usr/bin/env python3
"""Read-only complete C67 packet replay; no Magma job is launched.
Adapted for the repository: the checks of the packet manifest, of the .tex inputs and of the
Codex exit records are removed; the replays and the checks of the logs are unchanged.
"""
from pathlib import Path
import hashlib,json,subprocess,sys,time,re,ast
ROOT=Path(__file__).resolve().parents[1];C=ROOT/'common';P=ROOT/'P67';Q=P/'quotient45';D=P/'data';QD=D/'quotient45';LG=P/'logs'/'quotient45';t=time.monotonic()
log=(LG/'finite_validation_20261026/reemit_validated_data.log').read_text();assert 'NONIC67_STANDALONE_FINITE_DATA_PRODUCER_COMPLETE sources55=true targets77=true all_fibers_independently_complete=true'in log
assert len(re.findall(r'NONIC67_FRESH_TARGET_KEYS_COMPLETE fiber=\d degree=\d+ count11=true',log))==7
log=(LG/'modular_jbase_residue_check.log').read_text();assert 'NONIC67H3_NATIVE_MODULAR_NORM_CHECK_COMPLETE primes125 modulusdegree1500 heightbound748 all_exact_residues_agree=true'in log
log=(C/'logs'/'finite_remaining_census.log').read_text()
for d in[7,8,9]:
 line=next(x for x in log.splitlines()if f'N=67 degree={d} 'in x);assert'complete=true'in line
 orbits=ast.literal_eval(re.search(r'diamond_orbit_sizes=(\[.*?\])',line)[1]);js=ast.literal_eval(re.search(r'j_degrees=(\[.*?\])',line)[1]);assert sorted(zip(orbits,js))=={7:[(33,7)],8:[(33,8),(33,8)],9:[(11,3),(33,9),(33,9)]}[d]
log=(C/'logs'/'finite_remaining_primes_49.log').read_text();assert re.search(r'TORSION n=67 p=2.*= \[ 7, 8, 9 \]',log)
fi=json.loads((D/'FI_verified.json').read_text());assert fi['p']==67 and fi['H_index']==11 and fi['full_binary_rank']==22 and fi['prefix_configurations']==12233
r=json.loads((QD/'canonical_ring_sieve_all_verified.json').read_text());assert r['status']=='PASS'and r['cases']==550 and not r['unresolved']
assert {a['rank']for a in r['results']}=={132}and {a['positive_degree']for a in r['results']}=={56,64,72}
assert {d:sum(a['source_degree']==d for a in r['results'])for d in[7,8,9]}=={7:22,8:484,9:44}
commands=[[sys.executable,str(Q/p)]for p in['verify_q_identities.py','verify_finite_keys.py','verify_canonical_ring_tensor_readonly.py','verify_winding.py','verify_ranks_m4ri.py']]
commands += [[sys.executable,str(Q/'finite_validation_20261026/verify_reemitted_keys_readonly.py')],[sys.executable,str(Q/'canonical_ring_sieve.py'),'--readonly'],[sys.executable,str(Q/'verify_extended_q_and_sigma.py'),'--readonly'],[sys.executable,str(Q/'modular_jbase_norm.py'),'--readonly'],[sys.executable,str(C/'verify_FI_prefixes.py'),str(D/'FI_rows.txt')],[sys.executable,str(P/'certify_special9_readonly.py')]]
for cmd in commands:
 z=subprocess.run(cmd,capture_output=True,text=True,timeout=300);assert z.returncode==0,(cmd,z.stdout,z.stderr);print('CHECK_PASS',Path(cmd[1]).name,flush=True)
print('NONIC67_COMPLETE_READONLY_CERTIFICATE_PASS',json.dumps({'binary_matrices':550,'matrix_rank':132,'allbad_configurations':12233,'exact_CRT_base_prime_degree':1500,'claim':'No point of order67 over any degree9 number field','elapsed_seconds':time.monotonic()-t},sort_keys=True))
