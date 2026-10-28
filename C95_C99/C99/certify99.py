#!/usr/bin/env python3
"""Read-only complete degree9 C99 packet audit."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,math,re
import sympy as s
H=Path(__file__).resolve().parent;L=H.parent/'logs'/'C99';DD=H.parent/'data'/'C99';M=json.loads((DD/'manifest99.json').read_text())
assert M['degree']==9 and M['excluded_cyclic_orders']==[99] and M['complete']
log=(L/'winding99.log').read_text();assert log.rstrip().endswith('WINDING_CHARACTER_PROBE_DONE')
assert not re.search('Runtime error|User error|Assertion failed',log)
rr=[tuple(map(int,z)) for z in re.findall(r'WINDING N=99 order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+)',log)]
assert sum(int(s.totient(o))*d for o,c,d,w,cf in rr)==281
assert all(d==w for o,c,d,w,cf in rr if o>1)
assert pow(2,15,99) in [1,98] and all(pow(2,k,99) not in [1,98] for k in range(1,15))
assert 9%15 and 36<Fraction(43875,1024)
cm={int(n):int(d) for n,d,_ in re.findall(r'\[\*\s*(\d+),\s*(\d+),\s*(\d+),',(H.parent.parent/'common'/'data'/'dcm_list_all_min_orders_X1_10k.m').read_text())};assert cm[99]==60>9
F=json.loads((DD/'finite_trace_certificate.json').read_text());assert F['degree']==9 and F['N']==99 and len(F['rows'])==9
for r in F['rows']:
 e=r['e'];q=2**e;b=math.isqrt(4*q);ts=[q+1-h for h in range(99,q+b+2,99) if abs(h-q-1)<=b]
 ss={0,-2**((e+1)//2),2**((e+1)//2)} if e%2 else {0,-2**(e//2),2**(e//2),-2**(e//2+1),2**(e//2+1)}
 assert r['q']==q and r['candidate_traces']==ts and r['supersingular_traces']==sorted(ss)
 assert all(t%2==0 and t not in ss for t in ts) and r['admissible_traces']==[]
assert F['rows'][-1]['candidate_traces']==[18]
assert 'NONIC_C99_DIAMOND_PASS' in (L/'check_diamond99.log').read_text()
print('NONIC_C99_READONLY_CERTIFICATE_PASS degree=9 winding_genus=281 no_good_degrees=1..9 diamond_order=15 principal_bound=36 complete=true')
