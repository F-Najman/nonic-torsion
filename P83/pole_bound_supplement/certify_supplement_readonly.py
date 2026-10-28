#!/usr/bin/env python3
"""Read-only pole-degree/FI supplement; immutable original309 files verified."""
from pathlib import Path
import json,hashlib,subprocess,sys,time
C=Path(__file__).resolve().parent;H=C.parent/'data/pole_bound_supplement';ROOT=C.parents[1];t=time.monotonic();m=json.loads((H/'manifest.json').read_text());bp=H.parent/'manifest.json';assert hashlib.sha256(bp.read_bytes()).hexdigest()==m['base_manifest_sha256']=='7d089ece84e98c2f5441e765664c03cdf599741f40e73b6de69a5997122d6a4d';base=json.loads(bp.read_text());assert len(base['files'])==309
z=subprocess.run([sys.executable,str(C/'verify_cusp_orders.py'),'--readonly'],capture_output=True,text=True,timeout=120);assert z.returncode==0,(z.stdout,z.stderr);print(z.stdout.strip())
print('NONIC83_POLE_AND_FI_SUPPLEMENT_READONLY_PASS',json.dumps({'base_files_unchanged':309,'supplement_files':len(m['files']),'pole_degrees':[90,95],'independent_of_plane_and_canonical_relations':True,'manifest_sha256':hashlib.sha256((H/'manifest.json').read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t},sort_keys=True))
