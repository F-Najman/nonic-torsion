from pathlib import Path
import ast,json,sys
from fractions import Fraction
from sage.all import GF,PolynomialRing,EllipticCurve
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parents[1]/'data'/'C52';B=P.parent
def rd(p):
 tree=ast.parse(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'),mode='eval')
 def cv(n):
  if isinstance(n,ast.Constant) and isinstance(n.value,(int,bool)):return n.value
  if isinstance(n,ast.List):return [cv(x) for x in n.elts]
  if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -cv(n.operand)
  if isinstance(n,ast.BinOp) and isinstance(n.op,ast.Div):return Fraction(cv(n.left),cv(n.right))
  raise ValueError(ast.dump(n))
 return cv(tree.body)

a=rd(B/'N26_D9_threequotients_arithmetic.txt');records=a[3];orbits=a[8];assert [r[0] for r in a[4]]==[21,58]
row=records[57];assert row[0]==58;R=PolynomialRing(GF(3),'z');K=GF(3**9,'z',modulus=R(row[1]));z=K.gen()
def key(v):
 out=[]
 for i in range(9):out.append(tuple(int(x) for y in v for x in y._vector_()));v=[x**3 for x in v]
 return min(out)
example=rd(B/'explicit/example.txt');pol=R(example[0]);roots=pol.change_ring(K).roots(multiplicities=False);assert len(roots)==9
allkeys=[]
for root in roots:
 b=sum(K(t.numerator)*K(t.denominator)**(-1)*root**i for i,t in enumerate(example[2][0]));c=sum(K(t.numerator)*K(t.denominator)**(-1)*root**i for i,t in enumerate(example[2][1]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b);assert 26*Q==E(0) and 13*Q!=E(0) and 2*Q!=E(0) and T!=13*Q
 allkeys.append(key([b,c]));twos=[T,13*Q,T+13*Q];V=14*Q;ims=set()
 for ip in range(3):
  for iq in range(3):
   if ip==iq:continue
   for aa in range(1,7):
    U=twos[iq]+aa*V;PP=twos[ip];x=U[0]-PP[0];y=U[1];ims.add(key([x**3/y**2,(c+3*PP[0])*x*x/y**2]))
 assert ims==set(map(tuple,orbits[57]))
assert len(set(allkeys))==1 and len(set(map(tuple,orbits[57])))==12
cl=rd(B/'class_last2_arithmetic.txt');iv,rels,hgens,rows=cl[:4];assert iv==[2,14,70,1330,3990,11970,0]
H=hermite_normal_form(Matrix([v[:-1] for v in rels+hgens]).T);assert abs(int(H.det()))==350
flags=[]
for oid,pol,vals,vec,yes in rows:
 assert pol==records[oid-1][1] and vals==records[oid-1][2] and vec[-1]==0
 w=list(vec[:-1])
 for col in range(5,-1,-1):
  q=w[col]//int(H[col,col])
  for i in range(col+1):w[i]-=q*int(H[i,col])
 assert bool(yes)==(not any(w));flags.append([oid,bool(yes)])
assert flags==[[21,False],[58,True]] and cl[4]==[58]
assert sum(cl[5])==2*14*70*1330*3990*11970

out=dict(complete=True,degree=9,explicit_modular_degree=9,finite_surviving_orbit=58,explicit_embeddings_checked=9,automorphism_orbit_size=12,source_rational_image_index=350,last_two_class_memberships=flags)
if '--readonly' in sys.argv:assert json.loads((P/'uniqueness_sage_certificate.json').read_text())==out
else:(P/'uniqueness_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C52_EXPLICIT_C226_UNIQUE_ORBIT_LINK_SAGE_PASS',json.dumps(out),flush=True)
