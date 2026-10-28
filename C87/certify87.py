#!/usr/bin/env python3
"""Read-only complete degree9 C87 certificate; no Magma launch."""
from pathlib import Path
import json,hashlib,re,os,subprocess,sys
P=Path(__file__).resolve().parent;D=P/'data';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert (D/'finite_j_model_snapshot.m').read_bytes()==(D/'model_rebound_data.m').read_bytes()
lines=(D/'direct5_RR_cases.jsonseq').read_text().splitlines();rows=[json.loads(s)for s in lines]
assert len(rows)==64 and {(r['source_index'],r['target_index'])for r in rows}=={(i,j)for i in range(1,5)for j in range(1,17)} and all(r['dimension']==0 for r in rows)
assert 'Cases:=64;Survivors:=[];Complete:=true;' in(D/'direct5_RR_complete.m').read_text()
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
cmds=[[sys.executable,str(P/'direct5/certify_direct5.py')],[sys.executable,str(P/'modular_jbase_norm.py'),'--readonly'],[sys.executable,str(P/'canonical_replay/replay64.py'),'--readonly']]
for cmd in cmds:
 z=subprocess.run(cmd,env=env,capture_output=True,text=True);assert z.returncode==0,(cmd,z.stdout,z.stderr);print(z.stdout.strip(),flush=True)
print('NONIC87_COMPLETE_READONLY_CERTIFICATE_PASS degree9=true cases64=true native_RR_zero64=true independent_canonical_rank96=true no_classgroup_dependency=true',flush=True)
