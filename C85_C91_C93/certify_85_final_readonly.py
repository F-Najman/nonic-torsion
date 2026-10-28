#!/usr/bin/env python3
"""Portable module-path launcher for the frozen complete C85 replay."""
from pathlib import Path
import os,json,hashlib,subprocess,sys
C=Path(__file__).resolve().parent;R=C.parent;H=C/'C85';P=C/'data/C85/portability_supplement';m=json.loads((P/'launcher_manifest.json').read_text());assert m['portable_manifest_sha256']==hashlib.sha256((P/'manifest.json').read_bytes()).hexdigest()
env=os.environ.copy();env['PYTHONPATH']=str(H)+os.pathsep+env.get('PYTHONPATH','');env['PYTHONDONTWRITEBYTECODE']='1';z=subprocess.run([sys.executable,str(C/'certify_85_portable_readonly.py')],env=env);assert z.returncode==0
print('NONIC85_FINAL_READONLY_CERTIFICATE_PASS launcher_manifest_sha256='+hashlib.sha256((P/'launcher_manifest.json').read_bytes()).hexdigest())
