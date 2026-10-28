#!/usr/bin/env python3
"""Read-only full replay of the complete order73 exclusion."""
from pathlib import Path
import hashlib,json,sys,subprocess,time
ROOT=Path(__file__).resolve().parents[1];C=ROOT/'common';P=ROOT/'P73';Q=P/'quotient19';D=P/'data';LG=P/'logs';t0=time.monotonic()
log=(LG/'quotient19'/'native_control.log').read_text();assert 'NONIC73_RR_NATIVE_CONTROL_COMPLETE dimension=0'in log and 'NONIC73_SIGMA_NATIVE_CONTROL order4=true fixes_j=true'in log
commands=[[sys.executable,str(Q/s)]for s in ['verify_q_identities.py','verify_winding.py','verify_finite_keys.py','verify_canonical_ring_tensor_readonly.py','verify_multiplication_certificate_readonly.py','verify_ranks_m4ri.py']];commands+=[[sys.executable,str(C/'verify_FI_relations.py'),str(D/'FI_rows.txt')],[sys.executable,str(C/'replay_rational_two.py'),'73']]
for cmd in commands:
 r=subprocess.run(cmd,capture_output=True,text=True,timeout=180);assert r.returncode==0,(cmd,r.stdout,r.stderr);print('CHECK_PASS',Path(cmd[1]).name,flush=True)
print('NONIC73_COMPLETE_READONLY_CERTIFICATE_PASS',json.dumps({'binary_matrices':32,'matrix_shape':[90,164],'rank':90,'claim':'C73 excluded over every degree9 field','elapsed_seconds':time.monotonic()-t0},sort_keys=True))
