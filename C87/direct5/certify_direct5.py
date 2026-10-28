"""Read-only finite-source and Hecke-target replay. No global exclusion asserted."""
from pathlib import Path
import json,hashlib,os,resource,subprocess,sys
p=Path(__file__).resolve().parent
assert 'NONIC87_NATIVE5_TARGETS_COMPLETE' in (p.parent/'logs'/'direct5'/'isogeny5_targets.log').read_text()
def lim():
 resource.setrlimit(resource.RLIMIT_CPU,(600,600));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));os.nice(10)
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
r=subprocess.run([sys.executable,str(p/'replay5_sage.py'),'--check'],env=env,preexec_fn=lim,timeout=900);assert r.returncode==0
print('NONIC87_DIRECT5_READONLY_CERTIFICATE_PASS source_j_orbits=4 target_j_degrees=54 no_degree8_source=true')
