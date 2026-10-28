"""Read-only finite-fibre binding supplement; no Magma or Sage."""
from pathlib import Path
import json,hashlib,subprocess,sys
p=Path(__file__).resolve().parent;J=lambda q:json.loads(q.read_text());sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
h=p.parent/'data/C78_known39';assert sha(h/'manifest39.json')=='921fd3c267f053962ea1fb90c8936ac39be782114823b116c03c5c69e4e9ab05'
for cmd in[[sys.executable,str(p.parent/'C78_X139/census39_F512.py'),'--readonly'],[sys.executable,str(p/'bind117.py'),'--readonly']]:
 r=subprocess.run(cmd,capture_output=True,text=True);assert r.returncode==0,(r.stdout,r.stderr);print(r.stdout.strip())
print('C117_KNOWN39_FINITE_FIBRE_BINDING_CERTIFICATE_PASS')
