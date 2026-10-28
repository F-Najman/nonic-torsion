from math import gcd
from sympy import divisors,totient
import json
n=75;units=[u for u in range(1,n)if gcd(u,n)==1];rows=[]
for d in divisors(n):
 m=n//d;stab={u for u in units if u%m in({1%m,(-1)%m}if d<=2 else{1%m})};qdeg=len(units)//len(stab);f=1
 while pow(2,f,n)not in stab:f+=1
 rows.append({'d':int(d),'global_degree':qdeg,'finite_degree':f,'global_orbits':int(totient(d)*totient(m))//(2*qdeg)})
assert all(r['global_degree']==r['finite_degree']for r in rows if r['finite_degree']<=2)
assert all(7 in {u for u in units if u%(n//r['d']) in({1%(n//r['d']),(-1)%(n//r['d'])}if r['d']<=2 else{1%(n//r['d'])})}for r in rows if r['finite_degree']<=2)
print(json.dumps(rows));print('NONIC75_CUSP_COMPLEMENT_DEGREE2_PASS A7_kills_all=true')
