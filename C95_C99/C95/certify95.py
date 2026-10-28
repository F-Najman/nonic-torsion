#!/usr/bin/env python3
"""Read-only complete C95 exclusion audit, including the independent FI replay."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,math,re,runpy
import sympy as s
H=Path(__file__).resolve().parent;L=H.parent/'logs'/'C95';DD=H.parent/'data'/'C95';M=json.loads((DD/'manifest95.json').read_text())
assert M['degree']==9 and M['complete'] and M['excluded_cyclic_orders']==[95]
wind=(L/'winding95.log').read_text();assert wind.rstrip().endswith('WINDING_CHARACTER_PROBE_DONE')
rr=[tuple(map(int,z)) for z in re.findall(r'WINDING N=95 order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+)',wind)]
assert sum(int(s.totient(o))*d for o,c,d,w,cf in rr)==289 and [o for o,c,d,w,cf in rr if d!=w]==[12]
assert [(d,w) for o,c,d,w,cf in rr if o==1]==[(9,9)]
assert pow(11,3,95)==1 and all(pow(11,k,95) not in [1,94] for k in [1,2]) and Fraction(43875,1024)>36
log=(L/'check_diamond95.log').read_text();assert 'NONIC_95_FINITE_DIAMOND_IMAGE_KERNEL_PASS' in log
F=json.loads((DD/'finite_trace_certificate.json').read_text());assert F['N']==95 and F['degree']==9 and len(F['rows'])==9
for r in F['rows']:
 e=r['e'];q=2**e;b=math.isqrt(4*q);ts=[q+1-h for h in range(95,q+b+2,95) if abs(h-q-1)<=b]
 ss={0,-2**((e+1)//2),2**((e+1)//2)} if e%2 else {0,-2**(e//2),2**(e//2),-2**(e//2+1),2**(e//2+1)}
 assert r['q']==q and r['candidate_traces']==ts and r['supersingular_traces']==sorted(ss)
 assert all(t%2==0 and t not in ss for t in ts) and r['admissible_traces']==[]
assert F['rows'][7]['candidate_traces']==[-28] and F['rows'][8]['candidate_traces']==[38]
fi=(L/'x095_fi3.log').read_text();assert 'NONIC_95_X0_FI3_COMPLETE tested=20 exceptions=0' in fi
for fn in ['winding95.log','check_diamond95.log','x095_fi3.log']:
 assert not re.search('Runtime error|User error|Assertion failed',(L/fn).read_text())
cm={int(n):int(d) for n,d,_ in re.findall(r'\[\*\s*(\d+),\s*(\d+),\s*(\d+),',(H.parent.parent/'common'/'data'/'dcm_list_all_min_orders_X1_10k.m').read_text())};assert cm[95]==72>9
runpy.run_path(str(H/'replay_fi95.py'),run_name='__main__')
print('NONIC_C95_EXCLUSION_READONLY_PASS source_genus=289 X0_genus=9 all20_FI_ranks=3 complete=true')
