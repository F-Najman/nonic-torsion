#!/usr/bin/env python3
"""Build exact degree-nine global Hecke/gonality arithmetic and rank inputs."""
import ast,hashlib,json,math,re
from fractions import Fraction
from pathlib import Path
H=Path(__file__).resolve().parent;I=H/'data'/'inputs';L=H/'logs';ROOT=H.parent
cache=json.loads((ROOT/'common'/'data'/'lmfdb_newforms_cache.json').read_text())
old=json.loads((I/'degree8_large_certificate.json').read_text())
cm={int(n):int(d) for n,d,_ in re.findall(r'\[\*\s*(\d+),\s*(\d+),\s*(\d+),',(ROOT/'common'/'data'/'dcm_list_all_min_orders_X1_10k.m').read_text())}
assert sorted(cm)==list(range(1,10001))
def divisors(n):return [d for d in range(1,n+1) if n%d==0]
def phi(n):return sum(math.gcd(i,n)==1 for i in range(1,n+1))
def index(n):
 z=Fraction(n*n,2)
 for p in divisors(n):
  if p>1 and len(divisors(p))==2:z*=Fraction(p*p-1,p*p)
 assert z.denominator==1;return int(z)
def genus(n):
 z=1+Fraction(index(n),12)-Fraction(sum(phi(d)*phi(n//d) for d in divisors(n)),4)
 assert z.denominator==1;return int(z)
def order(n,a):
 assert math.gcd(n,a)==1
 b=a%n;o=1
 while b not in [1,n-1]:b=b*a%n;o+=1
 return o
def winding(n,path):
 s=path.read_text();assert s.rstrip().endswith('WINDING_CHARACTER_PROBE_DONE')
 assert not re.search(r'Runtime error|User error|Assertion failed|SystemError',s)
 rr=re.findall(r'WINDING N=(\d+) order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+) time=([0-9.]+)',s)
 rr=[dict(zip(['N','order','conductor','dimension','winding_dimension','coefficient_degree'],map(int,r[:6]))) for r in rr]
 assert rr and all(r['N']==n for r in rr)
 assert sum(phi(r['order'])*r['dimension'] for r in rr)==genus(n)
 m=re.search(r'WINDING_FINAL N=(\d+) genus=(\d+) total_dimension=(\d+) bad_orders=(\[.*?\]) admissible=(\[.*?\])',s)
 assert m and int(m[1])==n and int(m[2])==int(m[3])==genus(n)
 bad=ast.literal_eval(m[4]);assert bad==[r['order'] for r in rr if r['dimension']!=r['winding_dimension']]
 return rr,bad,ast.literal_eval(m[5])
rows=[];old_skipped=[]
for row in old['rows']:
 n=row['N'];a=row['a'];gon=Fraction(325*index(n),32768)
 if row['method']!='global3' or 9%(order(n,a))==0 or 9*row['kq']>=gon:
  old_skipped.append({'N':n,'reason':'local degree8 method or degree9 inequality/order not satisfied'});continue
 r={'N':n,'a':a,'diamond_order':order(n,a),'q':3,'kq':row['kq'],'index':index(n),'genus':genus(n),'gonality_lower_bound':str(gon),'degree9_principal_bound':9*row['kq'],'minimum_CM_degree':cm[n]}
 assert cm[n]>9 and n%3!=0 and r['diamond_order']==row['diamond_order']
 if row['rank_source']=='proved archived newforms':
  total=0;positive=[]
  for m in divisors(n):
   assert str(m) in cache
   for f in cache[str(m)]:
    assert f['analytic_rank_proved'] is True and f['analytic_rank'] is not None
    total+=len(divisors(n//m))*f['dim']
    if f['analytic_rank']>0:positive.append({'level':m,**f})
  assert total==genus(n)
  r.update(rank_source=row['rank_source'],newform_dimension_sum=total,positive_rank_forms=positive,
    nontrivial_positive_characters=sorted(set((f['level'],f['char_orbit_label']) for f in positive if f['char_order']>1)))
 else:
  rr,bad,admissible=winding(n,L/row['winding_log']);assert bad in [[],[1]] and a in admissible
  r.update(rank_source='degree-independent archived exact winding',winding_log=row['winding_log'],character_orbits=rr,bad_orders=bad)
 rows.append(r)
new_targets=[145,155,203,217,319,341,377,403,493,527,551,589,667,713,841,899,961]
unfinished=[]
for n in new_targets:
 p=L/f'winding{n}.log'
 fast=L/f'winding{n}_fast.log'
 if fast.exists() and fast.read_text().rstrip().endswith('WINDING_CHARACTER_PROBE_DONE'):p=fast
 if not p.exists() or not p.read_text().rstrip().endswith('WINDING_CHARACTER_PROBE_DONE'):
  unfinished.append(n);continue
 rr,bad,admissible=winding(n,p);gon=Fraction(325*index(n),32768)
 aa=[a for a in admissible if 9%order(n,a) and all(math.gcd(order(n,a),b)==1 for b in bad) and 9*(7 if a in [3,pow(3,-1,n)] else 8)<gon]
 if not aa:
  unfinished.append(n);continue
 a=3 if 3 in aa else min(aa);k=7 if a in [3,pow(3,-1,n)] else 8
 assert cm[n]>9 and n%3
 rows.append({'N':n,'a':a,'diamond_order':order(n,a),'q':3,'kq':k,'index':index(n),'genus':genus(n),'gonality_lower_bound':str(gon),'degree9_principal_bound':9*k,'minimum_CM_degree':cm[n],
 'rank_source':'fresh exact winding','winding_log':p.name,'character_orbits':rr,'bad_orders':bad})
rows.sort(key=lambda r:r['N'])
assert len(rows)==len(set(r['N'] for r in rows))
script=['SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);','AttachSpec("mdmagma/Magma/magma.spec");']
for r in rows:
 n=r['N'];a=r['a'];script.append(f'N:={n};a:={a};')
 if r['rank_source']=='proved archived newforms':
  for m,lab in r['nontrivial_positive_characters']:
   assert n==m
   script.extend([f'ii:=[j:j in [1..N-1]|GCD(j,N) eq 1 and ConreyCharacterOrbitLabel(N,j) eq "{m}.{lab}"];',
    'assert #ii gt 0 and &and[ConreyCharacterValue(N,j,a) eq 1:j in ii];',
    f'printf "NONIC_CHARACTER_KERNEL N=%o a=%o orbit={m}.{lab} count=%o exact=true\\n",N,a,#ii;'])
 else:
  # All characters whose order appears in the bad list are checked. This can
  # be stronger than necessary, but never assumes a good character was bad.
  bad=r['bad_orders'];script.append('badorders:={' + ','.join(map(str,bad)) + '};')
  script.append('chars:=[c:c in GaloisConjugacyRepresentatives(FullDirichletGroup(N))|IsEven(c) and Order(c) in badorders];')
  script.append('assert &and[Evaluate(c,a) eq 1:c in chars];')
  script.append('printf "NONIC_WINDING_BAD_CHARACTER_KERNEL N=%o a=%o badorders=%o checked_orbits=%o exact=true\\n",N,a,badorders,#chars;')
 script.append(f'printf "NONIC_GLOBAL_ROW N=%o a=%o diamond_order={r["diamond_order"]} d=9 principal_bound={r["degree9_principal_bound"]} exact=true\\n",N,a;')
script.extend(['print "NONIC_LARGE_CHARACTERS_COMPLETED";','quit;'])
(H/'check_nonic_characters.m').write_text('\n'.join(script)+'\n')
report={'degree':9,'method':'global Hecke/gonality plus free diamond action; degree-specific inequalities and rank certificates verified',
'excluded_cyclic_orders':[r['N'] for r in rows],'rows':rows,'unresolved_new_targets':unfinished,'old_rows_not_transferred':old_skipped}
(H/'data'/'large_global_certificate.json').write_text(json.dumps(report,indent=2)+'\n')
print('NONIC_LARGE_GLOBAL_ARITHMETIC_PASS exclusions='+str(len(rows))+' unfinished_new='+str(unfinished))
for r in rows:print(r['N'],r['a'],r['diamond_order'],r['degree9_principal_bound'],r['gonality_lower_bound'],r['minimum_CM_degree'])
