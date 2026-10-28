#!/usr/bin/env python3
"""Read-only complete degree9 C77 exclusion audit."""
from pathlib import Path
import hashlib,json,math,re,runpy
import sympy as s
H=Path(__file__).resolve().parent;
wind=(H/'logs'/'winding77.log').read_text();rr=[tuple(map(int,z)) for z in re.findall(r'NONIC77_WINDING order=(\d+) conductor=(\d+) relative_dim=(\d+) winding_dim=(\d+) coefficient_degree=(\d+)',wind)]
assert sum(int(s.totient(o))*d for o,c,d,w,cf in rr)==181
assert all(d==w for o,c,d,w,cf in rr if o>1)
assert [(d,w) for o,c,d,w,cf in rr if o==1]==[(7,6)]
assert 'NONIC77_POSITIVE_UPPER_T2_EXACTLY_ZERO [ 0 ]' in wind and 'NONIC77_FINITE_IMAGE_T2_AND_DIAMOND23_CERTIFIED' in wind
mod=(H/'logs'/'model_deck_rerun.log').read_text()
for marker in ['NONIC77_DIRECT_GENUS37_VERIFIED','NONIC77_GLOBAL_CUSPS_VERIFIED','NONIC77_DECK_DIFFERENTIALS_VERIFIED','NONIC77_BASIS_ALIGNMENT_PASS intertwiner_dimension=1 invertible=true','NONIC77_FI9_COMPLETE tested=3098 exceptions=0 complete=true']:assert marker in mod
geo=(H/'logs'/'validate_modular_cover.log').read_text()
m=re.search(r'NONIC77_COVER_IDENTIFICATION_DONE source_places=(\d+) remaining=1',geo);assert m and int(m[1])>=1
for marker in ['NONIC77_AS_CHARACTER_VERIFIED samples=5 all_Picard_mod2=true','NONIC77_FIXED_AS_CHARACTER_VALIDATED']:assert marker in geo
for text in [wind,mod,geo]:assert not re.search('Runtime error|User error|Assertion failed|FAILED',text)
F=json.loads((H/'data'/'finite_trace_certificate.json').read_text());assert F['N']==77 and F['degree']==9 and len(F['rows'])==9
for r in F['rows']:
 e=r['e'];q=2**e;b=math.isqrt(4*q);ts=[q+1-h for h in range(77,q+b+2,77) if abs(h-q-1)<=b]
 ss={0,-2**((e+1)//2),2**((e+1)//2)} if e%2 else {0,-2**(e//2),2**(e//2),-2**(e//2+1),2**(e//2+1)}
 assert r['q']==q and r['candidate_traces']==ts and r['supersingular_traces']==sorted(ss) and all(t%2==0 and t not in ss for t in ts)
assert [r['candidate_traces'] for r in F['rows'] if r['candidate_traces']]==[[-12],[26],[-26]]
runpy.run_path(str(H/'replay_fi9.py'),run_name='__main__')
print('NONIC_C77_EXCLUSION_READONLY_PASS source_genus=181 quotient_genus=37 all3098_FI_ranks=9 complete=true')
