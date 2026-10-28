#!/usr/bin/env python3
"""Read-only numerical/hash/complete-winding/modular-Krylov audit."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,re,math
import sympy as s
H=Path(__file__).resolve().parent.parent;L=H/'logs'/'addendum899';DD=H/'data'/'addendum899';M=json.loads((DD/'manifest.json').read_text());D=json.loads((DD/'certificate.json').read_text());n=D['N'];a=D['a'];order=D['diamond_order']
assert M['degree']==9 and M['complete'] and M['excluded_cyclic_orders']==[n]
log=(L/f'winding{n}_fast.log').read_text();assert log.rstrip().endswith('WINDING_CHARACTER_PROBE_DONE') and not re.search('Runtime error|User error|Assertion failed',log)
rr=[tuple(map(int,z)) for z in re.findall(r'WINDING N=\d+ order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+)',log)]
assert sum(int(s.totient(o))*d for o,c,d,w in rr)==D['genus']
assert [o for o,c,d,w in rr if d!=w]==D['bad_winding_orders']
assert pow(a,order,n) in [1,n-1] and all(pow(a,j,n) not in [1,n-1] for j in range(1,order))
assert 9%order and all(math.gcd(order,b)==1 for b in D['bad_winding_orders'])
idx=Fraction(n*n,2)
for p in s.factorint(n):idx*=Fraction(p*p-1,p*p)
g=1+idx/12-Fraction(sum(int(s.totient(d)*s.totient(n//d)) for d in s.divisors(n)),4)
assert idx==D['index'] and g==D['genus'] and D['degree9_principal_bound']==9*D['kq']<Fraction(325,32768)*idx==Fraction(D['gonality_lower_bound'])
cm={int(n):int(d) for n,d,_ in re.findall(r'\[\*\s*(\d+),\s*(\d+),\s*(\d+),',(H.parent/'common'/'data'/'dcm_list_all_min_orders_X1_10k.m').read_text())};assert cm[n]==D['minimum_CM_degree']>9
for fn in D['modular_Krylov_matrices']:
 d=json.loads((H/'data'/fn).read_text());p=d['p'];dim=d['d'];A=[d['A'][i*dim:(i+1)*dim] for i in range(dim)];v=d['v'];B=[]
 assert s.isprime(p)
 for i in range(dim):B.append(v);v=[sum(v[k]*A[k][j] for k in range(dim))%p for j in range(dim)]
 rank=0
 for j in range(dim):
  i=next((i for i in range(rank,dim) if B[i][j]%p),None)
  if i is None:continue
  B[i],B[rank]=B[rank],B[i];z=pow(B[rank][j],-1,p);B[rank]=[x*z%p for x in B[rank]]
  for k in range(rank+1,dim):
   z=B[k][j];B[k]=[(x-z*y)%p for x,y in zip(B[k],B[rank])]
  rank+=1
 assert rank==dim
assert f'NONIC_{n}_DIAMOND_KERNEL_PASS' in (L/'check.log').read_text()
print('NONIC_LARGEST_ADDENDUM_READONLY_PASS N='+str(n)+' genus='+str(D['genus'])+' matrices='+str(len(D['modular_Krylov_matrices']))+' complete=true')
