"""Exact stdlib-only binding of all order468 F512 level39 points to global witnesses."""
from pathlib import Path
import json,math,sys,hashlib,collections
p=Path(__file__).resolve().parents[1]/'data/C117_known39_binding';h=p.parent/'C78_known39';J=lambda q:json.loads(q.read_text());sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest();manifest=J(h/'manifest39.json')
assert sha(h/'manifest39.json')=='921fd3c267f053962ea1fb90c8936ac39be782114823b116c03c5c69e4e9ab05'
helper=J(h/'positive_divisors39.json');census=J(p.parent/'C78_X139/census39_F512.json');mod=(1<<9)+(1<<8)+1
assert helper['original_F512_polynomial']==census['field_modulus_coefficients']==[1,0,0,0,0,0,0,0,1,1]
def mul(a,b):
 c=0
 while b:
  if b&1:c^=a
  b>>=1;a<<=1
 while c.bit_length()>9:c^=mod<<(c.bit_length()-10)
 return c
def powr(a,n):
 c=1
 while n:
  if n&1:c=mul(c,a)
  a=mul(a,a);n>>=1
 return c
iv=[0]+[powr(a,510)for a in range(1,512)];assert all(mul(a,iv[a])==1 for a in range(1,512))
def trace(a):
 v=0
 for i in range(9):v^=a;a=mul(a,a)
 assert v in[0,1];return v
tr=[trace(a)for a in range(512)]
def jint(b,c):
 a1=1^c;b2=mul(a1,a1);b4=mul(a1,b);b6=mul(b,b);b8=mul(b,b6);c4=mul(b2,b2)
 disc=mul(mul(b2,b2),b8)^mul(b6,b6)^mul(mul(b2,b4),b6);assert disc
 return mul(mul(mul(c4,c4),c4),iv[disc])
def direct_count(b,c):
 a1=1^c;total=1
 for x in range(512):
  A=mul(a1,x)^b;rhs=mul(mul(x,x),x)^mul(b,mul(x,x))
  total+=1 if A==0 else(2 if tr[mul(rhs,iv[mul(A,A)])]==0 else 0)
 return total
bits=lambda v:sum(int(a)<<i for i,a in enumerate(v));geometric=set();levels=collections.Counter();orbits=[];direct=[]
for z in helper['known_positive_divisors']:
 dd,fp,bc=z['reduction_code'];assert dd==9 and fp==helper['original_F512_polynomial'];b,c=map(bits,bc);orbit=[]
 for i in range(9):
  pair=(powr(b,2**i),powr(c,2**i));orbit.append(pair)
 assert len(set(orbit))==9 and not(set(orbit)&geometric);geometric.update(orbit)
 for bb,cc in orbit:levels[jint(bb,cc)]+=1
 assert direct_count(b,c)==468;direct.append(468);orbits.append(orbit)
selected=census['selected_curves'];wanted={z['j']for z in selected if z['curve_order']==468};assert len(wanted)==3 and sorted(wanted)==[130,315,441]
assert set(levels)==wanted and all(levels[j]==12 for j in wanted)and len(geometric)==36
assert all(z['j_degree']==3 for z in selected if z['curve_order']==468);assert sum(z['curve_order']==546 for z in selected)==36
cp=[1]
for r in sorted(wanted):
 nc=[0]*(len(cp)+1)
 for i,a in enumerate(cp):nc[i]^=mul(a,r);nc[i+1]^=a
 cp=nc
assert len(cp)==4 and all(a in[0,1]for a in cp)
assert math.gcd(468,511)==1 and 2*3*13==78 and 9*13==117
bounds=[513-math.isqrt(2048),513+math.isqrt(2048)];assert bounds==[468,558]and[n for n in range(bounds[0],bounds[1]+1)if n%117==0]==[468]
out={'complete':True,'field_modulus_coefficients':helper['original_F512_polynomial'],'order468_j_invariants':sorted(wanted),'common_j_minimal_polynomial_over_F2':cp,'known_closed_nonic_points':4,'Frobenius_orbit_lengths':[9]*4,'distinct_geometric_marked_pairs':36,'marked_pairs_per_j':{str(j):levels[j]for j in sorted(wanted)},'direct_Tate_model_orders':direct,'finite_group_order':468,'gcd_order_with_q_minus_1':1,'order39_points_per_curve':24,'fine39_structures_per_curve':12,'all_order468_level39_points_exhausted':True,'Hasse_interval_F512':bounds,'only_multiple_of117_in_Hasse_interval':468,'scope':'Exact finite elliptic/level binding. Together with the separately proved X1(39) dyadic injectivity and lower-support exclusions this excludes C117; no C78 completeness claimed.'}
if'--readonly'in sys.argv:assert J(p/'binding117.json')==out
else:(p/'binding117.json').write_text(json.dumps(out,indent=2)+'\n')
print('C117_ALL_ORDER468_LEVEL39_BINDING_PASS known4_closed_points_geometric36_j130_315_441_each12=true')
