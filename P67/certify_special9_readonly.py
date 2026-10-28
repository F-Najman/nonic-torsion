#!/usr/bin/env python3
"""Read-only replay of the exceptional C67 branch, not a full exclusion."""
from pathlib import Path
import hashlib,json,subprocess,sys
P=Path(__file__).resolve().parent
r=subprocess.run([sys.executable,str(P/'replay_special9.py')],capture_output=True,text=True,timeout=30);assert r.returncode==0,(r.stdout,r.stderr)
print('NONIC67_SPECIAL9_PARTIAL_CERTIFICATE_PASS',json.dumps({'fixed_special_points':11,'primary_cubic_point_theorem':'BGG arXiv2102.08236v2 Theorem1.1 and Section5.2.5','both_cubic_fields_totally_ramified_at_2':True,'complete_C67_exclusion':False},sort_keys=True))
