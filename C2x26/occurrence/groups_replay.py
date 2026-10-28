"""Independent primitive levels and exhaustive F49/F343 point counts."""
from pathlib import Path
import ast,sys,importlib.util
P=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/'data'/'occurrence'
BASE=Path(__file__).resolve().parents[2]

def module(name,path):
 sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
m3=module('c226_field3',BASE/'common/finite_arithmetic_char3.py')
m7=module('c226_field7',BASE/'C2x30/finite_arithmetic.py')
s=(P/'finite_geometric_arithmetic.txt').read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False');labels,av,am,bv,cv,cm,records=ast.literal_eval(s)
assert len(records)==2 and [r[0] for r in records]==[3,7]
counts=[]
for row in records:
 p,da,dr,dl,ds,ms,localdata=row[:7];assert da==2 and dr==1 and dl==1
 assert ds==[d[0] for d in localdata] and ms==[d[1] for d in localdata] and all(m==1 for m in ms)
 assert sum(d*m for d,m in zip(ds,ms))==9
 for de,m,n,pol,vals in localdata:
  mm=m3 if p==3 else m7;K=mm.Field(pol);assert K.n==de
  b,c=[mm.evaluate(v,K.z) for v in vals];E=mm.Curve(K,b,c);P2=(K.zero,K.zero);Q=(b,b)
  E.check(P2);E.check(Q);assert E.mul(2,P2) is None and E.mul(26,Q) is None and E.mul(13,Q) is not None and E.mul(2,Q) is not None and P2!=E.mul(13,Q)
  assert min(d for d in range(1,de+1) if de%d==0 and b**(p**d)==b and c**(p**d)==c)==de
  dd=b*(1-b-c);assert dd!=0 and c*c-4*dd!=0
  jj=256*(c*c-3*dd)**3/(dd*dd*(c*c-4*dd));jd=min(d for d in range(1,de+1) if de%d==0 and jj**(p**d)==jj)
  if p==3:
   assert de==9 and n==19604 and jd==3
  else:
   total=1
   for k in range(K.q):
    digits=[(k//(7**i))%7 for i in range(de)];x=mm.Elt(K,digits);y2=x*x*x+E.a2*x*x+E.a4*x
    if y2==0:total+=1
    elif y2**((K.q-1)//2)==1:total+=2
   assert total==n
   counts.append((de,n));print('C226_INDEPENDENT_FINITE_GROUP',p,de,n,'j_degree',jd)
assert sorted(counts)==[(2,52),(2,52),(2,52),(3,364)]
assert records[0][4:6]==[[9],[1]] and sorted(records[1][4])==[2,2,2,3]
print('C226_INDEPENDENT_TORSION_BOUND_PASS full_group_order_bound=52 lower_subgroup=2x26 unramified_good_prime=7 residue_degree=2')
print('C226_INDEPENDENT_LEVEL_AND_DEGREE_PASS source_closed_degree=9 j_degree=3')
