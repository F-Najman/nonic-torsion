#!/usr/bin/env python3
"""Read-only complete numerical, hash and modular-matrix replay for C841."""
from pathlib import Path
import hashlib,json,re,math
from fractions import Fraction
import sympy as s
H=Path(__file__).resolve().parent.parent;L=H/'logs'/'addendum841';DD=H/'data'/'addendum841';M=json.loads((DD/'manifest841.json').read_text());D=json.loads((DD/'certificate841.json').read_text())
log=(L/'winding841_fast.log').read_text();assert log.rstrip().endswith('WINDING_CHARACTER_PROBE_DONE')
assert not re.search('Runtime error|User error|Assertion failed',log)
rr=[tuple(map(int,z)) for z in re.findall(r'WINDING N=841 order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+)',log)]
assert sum(int(s.totient(o))*d for o,c,d,w in rr)==28834
assert [o for o,c,d,w in rr if d!=w]==[1,7]
assert pow(12,58,841) in [1,840] and all(pow(12,j,841) not in [1,840] for j in range(1,58))
assert math.gcd(58,7)==1 and 9%58 and 72<Fraction(325*353220,32768)==Fraction(D['gonality_lower_bound'])
cm={int(n):int(d) for n,d,_ in re.findall(r'\[\*\s*(\d+),\s*(\d+),\s*(\d+),',(H.parent/'common'/'data'/'dcm_list_all_min_orders_X1_10k.m').read_text())};assert cm[841]==406>9
for fn in D['modular_Krylov_matrices']:
 d=json.loads((H/'data'/fn).read_text());p=d['p'];n=d['d'];A=[d['A'][i*n:(i+1)*n] for i in range(n)];v=d['v'];B=[]
 assert s.isprime(p)
 for i in range(n):B.append(v);v=[sum(v[k]*A[k][j] for k in range(n))%p for j in range(n)]
 rank=0
 for j in range(n):
  i=next((i for i in range(rank,n) if B[i][j]%p),None)
  if i is None:continue
  B[i],B[rank]=B[rank],B[i];z=pow(B[rank][j],-1,p);B[rank]=[x*z%p for x in B[rank]]
  for k in range(rank+1,n):
   z=B[k][j];B[k]=[(x-z*y)%p for x,y in zip(B[k],B[rank])]
  rank+=1
 assert rank==n and n in [71,72]
assert 'NONIC_841_DIAMOND_KERNEL_PASS order=58' in (L/'check841.log').read_text()
print('NONIC_C841_READONLY_PASS genus=28834 full_Krylov_ranks=72,71 diamond_order=58 complete=true')
