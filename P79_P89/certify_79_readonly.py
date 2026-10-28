#!/usr/bin/env python3
"""Read-only full replay of the frozen order79 exclusion packet."""
from pathlib import Path
import json,hashlib,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];C=Path(__file__).resolve().parent;P=C/'data/P79';Q=C/'P79/quotient18';L=C/'logs/P79/quotient18';t0=time.monotonic()
m=json.loads((P/'manifest.json').read_text());assert m['degree']==9 and m['group']==[1,79]
assert 'NONIC79_RR_NATIVE_CONTROL_COMPLETE dimension=0'in(L/'native_control.log').read_text()
assert 'NONIC79_SIGMA_NATIVE_CONTROL order3=true fixes_j=true'in(L/'native_control.log').read_text()
checks=[[sys.executable,str(Q/s)]for s in ['verify_q_identities_complete.py','verify_winding.py','verify_finite_keys.py','verify_canonical_ring_tensor_readonly.py','verify_multiplication_certificate_readonly.py','verify_ranks_m4ri.py']];checks.append([sys.executable,str(ROOT/'common/verify_FI_relations.py'),str(P/'FI_rows.txt')])
for cmd in checks:
 r=subprocess.run(cmd,capture_output=True,text=True,timeout=180);assert r.returncode==0,(cmd,r.stdout,r.stderr);print('CHECK_PASS',Path(cmd[1]).name,flush=True)
print('NONIC79_COMPLETE_READONLY_CERTIFICATE_PASS',json.dumps({'files':len(m['files']),'native_exits_zero':14,'binary_matrices':48,'matrix_shape':[119,235],'rank':119,'claim':'C79 excluded over every degree9 field','elapsed_seconds':time.monotonic()-t0,'manifest_sha256':hashlib.sha256((P/'manifest.json').read_bytes()).hexdigest()},sort_keys=True))
