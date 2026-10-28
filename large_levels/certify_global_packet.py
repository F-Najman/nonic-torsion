#!/usr/bin/env python3
"""Read-only hash, rank-input, inequality and modular-Krylov proof audit."""
import ast,hashlib,json,math,re
from fractions import Fraction
from pathlib import Path
import sympy as s
H=Path(__file__).resolve().parent;L=H/'logs';DD=H/'data';ROOT=H.parent
M=json.loads((DD/'global_manifest.json').read_text());D=json.loads((DD/'large_global_certificate.json').read_text())
assert M['degree']==D['degree']==9 and M['complete']
assert M['excluded_cyclic_orders']==D['excluded_cyclic_orders']
charlog=(L/'check_nonic_characters.log').read_text()
assert charlog.rstrip().endswith('NONIC_LARGE_CHARACTERS_COMPLETED') and not re.search('Runtime error|Assertion failed|User error',charlog)
cache=json.loads((ROOT/'common/data/lmfdb_newforms_cache.json').read_text())
cm={int(n):int(d) for n,d,_ in re.findall(r'\[\*\s*(\d+),\s*(\d+),\s*(\d+),',(ROOT/'common/data/dcm_list_all_min_orders_X1_10k.m').read_text())}
matrices=set()
for r in D['rows']:
 n,a=r['N'],r['a'];assert n%3 and math.gcd(n,a)==1
 idx=Fraction(n*n,2)
 for p in s.factorint(n):idx*=Fraction(p*p-1,p*p)
 g=1+idx/12-Fraction(sum(int(s.totient(d)*s.totient(n//d)) for d in s.divisors(n)),4)
 assert idx==r['index'] and g==r['genus']
 o=int(s.n_order(a,n))
 if o%2==0 and pow(a,o//2,n)==n-1:o//=2
 assert o==r['diamond_order'] and 9%o!=0
 k=7 if a in [3,pow(3,-1,n)] else 8
 assert k==r['kq'] and 9*k==r['degree9_principal_bound']<Fraction(325,32768)*idx
 assert str(Fraction(325,32768)*idx)==r['gonality_lower_bound'] and cm[n]==r['minimum_CM_degree']>9
 assert f'NONIC_GLOBAL_ROW N={n} a={a} diamond_order={o} d=9 principal_bound={9*k} exact=true' in charlog
 if r['rank_source']=='proved archived newforms':
  total=0;pos=[]
  for lev in s.divisors(n):
   for f in cache[str(lev)]:
    assert f['analytic_rank_proved'] and f['analytic_rank'] is not None
    total+=len(s.divisors(n//lev))*f['dim']
    if f['analytic_rank']>0:pos.append({'level':int(lev),**f})
  assert total==g and pos==r['positive_rank_forms']
 else:
  log=(L/r['winding_log']).read_text();assert log.rstrip().endswith('WINDING_CHARACTER_PROBE_DONE')
  assert not re.search('Runtime error|Assertion failed|User error|SystemError',log)
  rr=re.findall(r'WINDING N=(\d+) order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+)',log)
  assert sum(int(s.totient(int(z[1])))*int(z[3]) for z in rr)==g
  assert all(int(z[0])==n for z in rr)
  bad=[int(z[1]) for z in rr if int(z[3])!=int(z[4])]
  assert bad==r['bad_orders']
  if r['rank_source']=='fresh exact winding':assert all(math.gcd(o,b)==1 for b in bad)
  for name in re.findall(r'WINDING_KRYLOV_MATRIX N=\d+ chi_index=\d+ order=\d+ path=(\S+)',log):
   matrices.add(Path(name).name)
for name in sorted(matrices):
 d=json.loads((DD/name).read_text());p=d['p'];dim=d['d'];v=d['v'];a=d['A']
 assert s.isprime(p) and len(v)==dim and len(a)==dim*dim
 A=[a[i*dim:(i+1)*dim] for i in range(dim)];B=[];u=v[:]
 for _ in range(dim):
  B.append(u);u=[sum(u[i]*A[i][j] for i in range(dim))%p for j in range(dim)]
 rank=0
 for j in range(dim):
  ii=next((i for i in range(rank,dim) if B[i][j]%p),None)
  if ii is None:continue
  B[ii],B[rank]=B[rank],B[ii];inv=pow(B[rank][j],-1,p);B[rank]=[(z*inv)%p for z in B[rank]]
  for i in range(rank+1,dim):
   z=B[i][j];B[i]=[(x-z*y)%p for x,y in zip(B[i],B[rank])]
  rank+=1
 assert rank==dim,name
print('NONIC_GLOBAL_PACKET_READONLY_PASS exclusions='+str(len(D['rows']))+' modular_Krylov_replays='+str(len(matrices))+' complete=true')
