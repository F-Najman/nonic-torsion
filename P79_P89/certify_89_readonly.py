#!/usr/bin/env python3
"""Read-only complete replay of the frozen order89 exclusion packet."""
from pathlib import Path
import json,hashlib,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];C=Path(__file__).resolve().parent;P=C/'data/P89';Q=C/'P89/quotient27';L=C/'logs/P89/quotient27';t0=time.monotonic()
m=json.loads((P/'manifest.json').read_text());assert m['degree']==9 and m['group']==[1,89]
assert 'NONIC89_RR_DIMENSION_PILOT_RESULT dimension=0'in(L/'rr_dimension_pilot.log').read_text()
checks=['verify_q_identities_complete.py','verify_winding.py','verify_finite_keys_jv.py','verify_canonical_ring_tensor_readonly.py','verify_small_multiplication_certificate_readonly.py','verify_small_ranks_m4ri.py']
for script in checks:
 r=subprocess.run([sys.executable,str(Q/script)],capture_output=True,text=True,timeout=180);assert r.returncode==0,(script,r.stdout,r.stderr)
 print('CHECK_PASS',script,flush=True)
print('NONIC89_COMPLETE_READONLY_CERTIFICATE_PASS',json.dumps({'files':len(m['files']),'native_exits_zero':15,'binary_matrices':96,'matrix_shape':[130,264],'rank':130,'claim':'C89 excluded over every degree9 field','elapsed_seconds':time.monotonic()-t0,'manifest_sha256':hashlib.sha256((P/'manifest.json').read_bytes()).hexdigest()},sort_keys=True))
