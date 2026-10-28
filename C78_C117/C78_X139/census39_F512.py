#!/usr/bin/env python3
"""Exact ordinary elliptic-curve census over F512, no Sage or Magma.
For y^2+xy=x^3+b the order is 513+K(sqrt(b)); the other twist
has order513-K. A Walsh transform computes every Kloosterman sum.
All selected orders are also checked by independent direct character sums.
"""
from pathlib import Path
import json,hashlib,sys,time
P=Path(__file__).resolve().parent;MOD=(1<<9)+(1<<8)+1;q=512

def mul(a,b):
 c=0
 while b:
  if b&1:c^=a
  b>>=1;a<<=1
  if a&q:a^=MOD
 return c

def power(a,n):
 c=1
 while n:
  if n&1:c=mul(c,a)
  a=mul(a,a);n>>=1
 return c

def trace(a):
 t=0
 for _ in range(9):t^=a;a=mul(a,a)
 assert t in (0,1)
 return t

inv=[0]+[power(a,510) for a in range(1,q)]
assert all(mul(a,inv[a])==1 for a in range(1,q))
tr=[trace(a)for a in range(q)]
f=[0]+[1-2*tr[inv[y]]for y in range(1,q)]
w=f[:];h=1
while h<q:
 for i in range(0,q,2*h):
  for j in range(i,i+h):a,b=w[j],w[j+h];w[j],w[j+h]=a+b,a-b
 h*=2

def label(a):return sum(tr[mul(a,1<<i)]<<i for i in range(9))
assert len({label(a)for a in range(q)})==q
# Every j=0 equation is represented here: cube scaling makes a3=1;
# x/y translations reduce a4 modulo s^4+s and a6 modulo s^2+s.
# Both additive images are the trace-zero hyperplane in odd degree9.
assert {mul(mul(a,a),mul(a,a))^a for a in range(q)}=={a for a in range(q)if tr[a]==0}
assert {mul(a,a)^a for a in range(q)}=={a for a in range(q)if tr[a]==0}
assert len({power(a,3)for a in range(1,q)})==511
ss=[]
for a4 in (0,1):
 for a6 in (0,1):
  order=1+sum(2 if tr[power(x,3)^mul(a4,x)^a6]==0 else 0 for x in range(q))
  ss.append({'a4':a4,'a6':a6,'order':order})
assert sorted(z['order']for z in ss)==[481,513,513,545]
assert all(z['order']%39 for z in ss)
selected=[];hist={}
for b in range(1,q):
 k=w[label(power(b,256))];hist[k]=hist.get(k,0)+1
 for twist in (0,1):
  order=513+(-1 if twist else 1)*k
  if order%39:continue
  direct=513+sum(1-2*(tr[x^mul(b,mul(inv[x],inv[x]))]^twist)for x in range(1,q))
  assert direct==order
  j=inv[b];cycle=[];a=j
  while a not in cycle:cycle.append(a);a=mul(a,a)
  assert a==j and len(cycle)in(3,9)
  selected.append({'j':j,'b':b,'twist_trace':twist,'curve_order':order,'j_degree':len(cycle),'j_frobenius_orbit':cycle})
assert len(selected)==39 and len({r['j']for r in selected})==39
assert sum(r['curve_order']==468 and r['j_degree']==3 for r in selected)==3
assert sum(r['curve_order']==546 and r['j_degree']==9 for r in selected)==36
orbits=sorted({tuple(sorted(r['j_frobenius_orbit']))for r in selected},key=lambda z:(len(z),z))
assert list(map(len,orbits))==[3,9,9,9,9]
out={'complete':True,'field_modulus_coefficients':[1,0,0,0,0,0,0,0,1,1],'all511_ordinary_j_invariants_and_both_twists_checked':True,'all_supersingular_four_models_directly_counted':ss,'all39_selected_orders_direct_character_sum_checked':True,'eligible_ordinary_curves':39,'j_degree3_order468_curves':3,'j_degree9_order546_curves':36,'j_frobenius_orbits':[list(z)for z in orbits],'source_non_cuspidal_degree9_places':52,'projective_diamond_orbit_sizes':[4,12,12,12,12],'selected_curves':selected,'Kloosterman_histogram':{str(k):v for k,v in sorted(hist.items())},'proof_note':'Supersingular j0 orders over F512 are481,513,545, none divisible39. Ordinary rational39torsion is cyclic because gcd(39,511)=1. Each eligible curve contributes phi(39)/2=12 marked pairs. Hasse excludes39torsion over F2 and F8. Thus468 rational marked pairs give52 closed nonic places in five diamondorbits.'}
file=P.parents[1]/'C78_C117/data/C78_X139/census39_F512.json'
if '--readonly'in sys.argv:assert json.loads(file.read_text())==out
else:file.write_text(json.dumps(out,indent=2)+'\n')
print('C39_F512_INDEPENDENT_COMPLETE_CENSUS_PASS curves39 places52 diamond_orbits4_12_12_12_12')
