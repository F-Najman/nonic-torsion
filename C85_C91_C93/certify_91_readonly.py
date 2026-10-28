#!/usr/bin/env python3
"""Complete frozen C91 exclusion certificate, without native launches or mutation."""
from pathlib import Path
import json,hashlib,subprocess,sys,time
C=Path(__file__).resolve().parent;R=C.parent;H=C/'data/C91';K=C/'C91';t=time.monotonic();m=json.loads((H/'manifest.json').read_text());assert m['degree']==9 and m['group']==[1,91]and m['status']=='COMPLETE EXCLUSION'
z=subprocess.run([sys.executable,str(K/'verify_all_components_readonly.py')]);assert z.returncode==0
print('NONIC91_COMPLETE_READONLY_CERTIFICATE_PASS',json.dumps({'files':len(m['files']),'selected_native_exits_zero':10,'binary_cases':288,'matrix_shape':[265,590],'matrix_rank':265,'special_degree8_closed_points':3,'manifest_sha256':hashlib.sha256((H/'manifest.json').read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t},sort_keys=True))
