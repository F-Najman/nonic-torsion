"""Read-only known nonic39 global/full-torsion and positive-reduction binding audit."""
from pathlib import Path
import json,hashlib,shutil,subprocess,sys
p=Path(__file__).resolve().parent;d=p.parents[1]/'C78_C117/data/C78_known39';J=lambda f:json.loads((d/f).read_text());sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest();m=J('manifest39.json')
b=J('positive_divisors39.json');assert b['global_closed_nonic_divisors']==4 and len(b['known_positive_divisors'])==4
assert J('independent39.json')['full_global_torsion_order']==39
r=subprocess.run([sys.executable,str(p/'replay39.py'),'--check'],capture_output=True,text=True);assert r.returncode==0,(r.stdout,r.stderr);print(r.stdout.strip())
print('NONIC78_KNOWN39_READONLY_CERTIFICATE_PASS',len(m['files']),'files; known4nonic39divisors only, no C78exclusion asserted')
