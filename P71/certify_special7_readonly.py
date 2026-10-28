#!/usr/bin/env python3
"""Read-only replay of a partial C71 branch exclusion. No Magma is launched."""
from pathlib import Path
import hashlib,json,subprocess,sys
P=Path(__file__).resolve().parent;LG=P/'logs'
assert 'degree7_closed=40 fixed_by_delta32=5 fixed_j_degree1=true diamond_orbits=2 complete=true' in (LG/'special7_finite_deck.log').read_text()
assert 'RATIONAL_TWO_DONE N=71 kernel_dimension=0 diamond_trivial=true' in (LG/'special7_rational_two.log').read_text()
r=subprocess.run([sys.executable,str(P/'replay_special7_two.py')],capture_output=True,text=True,timeout=30);assert r.returncode==0,(r.stdout,r.stderr)
print('NONIC71_SPECIAL7_PARTIAL_CERTIFICATE_PASS',json.dumps({'integral_homology_rank':352,'rational_two_torsion_upper_bound_dimension':0,'remaining_branches':['generic NC7 plus rational cusps of degree2','NC8 plus rational cusp'],'complete_C71_exclusion':False},sort_keys=True))
