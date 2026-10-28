"""Read-only exact Sage replay and manifest verifier."""
from pathlib import Path
import os,sys,json,hashlib,subprocess,resource
p=Path(__file__).resolve().parent
def limits():
 resource.setrlimit(resource.RLIMIT_CPU,(300,300));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));os.nice(10)
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
r=subprocess.run([sys.executable,str(p/'verify_known31.py'),'--check'],env=env,preexec_fn=limits,timeout=600)
assert r.returncode==0,r.returncode
print('KNOWN31_READONLY_CERTIFICATE_PASS full_torsion=31 all_five_survivors_bound=true')
