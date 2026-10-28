from pathlib import Path
import ast,json,sys
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
BASE=Path(__file__).resolve().parents[1]/'data';P=BASE/'diagonal_quotient'

def rd(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
class Quotient:
 def __init__(self,R,hs):
  self.n=len(R[0])-1;self.H=hermite_normal_form(Matrix([r[:-1] for r in R+hs]).T);self.index=abs(int(self.H.det()))
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0;w=list(v[:-1])
  for j in range(self.n-1,-1,-1):
   q=w[j]//int(self.H[j,j])
   for i in range(j+1):w[i]-=q*int(self.H[i,j])
  return tuple(w)
def fkey(vs,D):
 out=[]
 for i in range(D):out.append(tuple(int(x) for v in vs for x in v._vector_()));vs=[v**5 for v in vs]
 return min(out)
hom,images,tests,remaining=rd(P/'quotient_filter_arithmetic.txt');iv,rels,hc,acs,a,aa,ker,kiv,iiv,hq,hqiv=hom
assert iv==[3,3,3,24,24,192,0] and kiv==[3,3,3,24,24,24] and iiv==[8] and all(not any(v) for v in acs)
G=Quotient(rels,[]);H=Quotient(rels,hc);assert G.index==2985984 and H.index==8
# Native ker vectors use the kernel subgroup basis; the graph proves its ambient image equals H.
assert len({H.key([n*v for v in a]) for n in range(8)})==8
assert G.key([8*v for v in aa])==(0,)*6 and all(G.key([n*v for v in aa])!=(0,)*6 for n in [1,2,4])
# Independent graph-lattice proof that cusp images0 plus a->aa define a hom.
r=G.n;rr=[v[:-1] for v in rels];gens=hc+[a];outs=acs+[aa]
rows=[v+[0]*r for v in rr]+[[0]*r+v for v in rr]+[u[:-1]+v[:-1] for u,v in zip(gens,outs)]
graph=hermite_normal_form(Matrix(rows).T);assert abs(int(graph.det()))==G.index
assert Quotient(rels,gens).index==1 and Quotient(rels,outs).index==G.index//8
first=rd(BASE/'filter_quotient18_linear_arithmetic.txt');single,small,survivors=first[2],first[3],first[6];assert len(survivors)==len(tests)==2064
own={(v[0],v[1]):v for v in single};sm={(v[0],v[1]):v for v in small};R=PolynomialRing(GF(5),'z');fields={d:GF(5**d,'z',modulus=R(pol)) for d,pol in first[1]};cache={}
for row in images:
 d,key,e,w,cl,projkey=row;F=fields[d];z=F.gen();b=sum(F(v)*z**j for j,v in enumerate(key[:d]));c=sum(F(v)*z**j for j,v in enumerate(key[d:]));E=EllipticCurve(F,[0,c,0,b*(1-b-c),0]);Q=E(b,b);Q12=3*Q;T=12*Q;R9=4*Q
 lam=Q12[1]/Q12[0];xp=R9[0]+(R9+T)[0]+(R9-T)[0]-2*T[0];vals=[Q12[0]/lam**2,c/lam**2,xp/lam**2]
 assert fkey(vals,d)==tuple(projkey);de=min(int(j) for j in divisors(d) if all(v**(5**int(j))==v for v in vals));assert e==de and w==d//e and w in [1,3]
 cache[tuple([d]+key)]=cl
expected=[]
for k,desc,cl,yes in tests:
 assert desc==survivors[k-1];rows=[own[desc[1],desc[2]]]
 if desc[0]>=2:rows.append(sm[desc[3],desc[4]])
 if desc[0]==3:rows.append(sm[desc[5],desc[6]])
 summed=[0]*7
 for row in rows:summed=[a+b for a,b in zip(summed,cache[tuple([row[0]]+row[2])])]
 assert G.key(summed)==G.key(cl) and bool(yes)==(H.key(cl)==(0,)*6)
 if yes:expected.append(desc)
assert expected==remaining
out=dict(complete=True,degree=9,prime=5,quotient_genus=9,finite_picard_order=G.index,finite_A7_image_order=8,finite_A7_kernel_invariants=kiv,rational_image_upper_index=8,normalized_quotient_images_checked=len(images),supports=2064,remaining=len(remaining))
if '--readonly' in sys.argv:assert json.loads((P/'quotient_sage_certificate.json').read_text())==out
else:(P/'quotient_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C236_SAGE_QUOTIENT_HECKE_KERNEL_AND_ALL_PROJECTIONS_PASS',json.dumps(out),flush=True)
