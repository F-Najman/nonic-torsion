from pathlib import Path
import ast,json,sys
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent/'data'
def data(n):return ast.literal_eval((P/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
class Quotient:
 def __init__(self,R,H):self.n=len(R[0])-1;self.M=hermite_normal_form(Matrix([r[:-1] for r in R+H]).T);self.index=abs(int(self.M.det()))
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0;w=list(v[:-1])
  for j in range(self.n-1,-1,-1):
   q=w[j]//int(self.M[j,j])
   for i in range(j+1):w[i]-=q*int(self.M[i,j])
  return tuple(w)
def norm(E,T,Q):x=Q[0]-T[0];y=Q[1];return [x**3/y**2,(E.a2()+3*T[0])*x*x/y**2]
def iso3(E,T,U):
 def point(Q):
  S=Q+T;V=Q-T
  return (Q[0]+S[0]+V[0]-2*T[0],Q[1]+S[1]+V[1])
 pts=[point(Q) for Q in U];(x0,y0),(x1,y1)=pts
 a2=E.a2();v0=y0*y0-x0**3-a2*x0*x0;v1=y1*y1-x1**3-a2*x1*x1;a4=(v1-v0)/(x1-x0);a6=v0-a4*x0;EP=EllipticCurve(E.base_field(),[0,a2,0,a4,a6]);out=[EP(x,y) for x,y in pts]
 assert 2*out[1]==EP(point(2*U[1]));return EP,out
f=data('filter_quotient14_linear_arithmetic.txt');g,fs,single,small,pairs,triples,survivors,coverage=f;iv,rels,cd,cc,hg,hinv=g;assert iv==[2,2,12,36,0] and hinv==[2,2,6,18] and triples==[];H=Quotient(rels,hg);G=Quotient(rels,[]);assert H.index==4 and all(H.key(v)==(0,0,0,0) for v in cc)
R=PolynomialRing(GF(5),'z');fields={d:GF(5**d,'z',modulus=R(pol)) for d,pol in fs};own={};sm={}
def source(row):
 d,i,ky=row[:3];K=fields[d];z=K.gen();b=sum(K(v)*z**j for j,v in enumerate(ky[:d]));c=sum(K(v)*z**j for j,v in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);return E,E(0,0),E(b,b)
def degree(pair,d):return min(int(e) for e in divisors(d) if all(v**(5**int(e))==v for v in pair))
for is_small,rows in [(False,single),(True,small)]:
 for row in rows:
  d,i,ky,e,w,cl=row[:6];E,T,Q=source(row);de=degree(norm(E,T,3*Q),d);assert e==de and w==d//e;(sm if is_small else own)[d,i]=row
  if not is_small:assert bool(row[6])==(H.key(cl)==(0,0,0,0))
# Link every single row to the independently audited complete source list.
for d in range(4,10):
 raw=data('M2N42_D9_direct_full%d_arithmetic.txt'%d);assert [r[2] for r in single if r[0]==d]==[ky for de,ky in raw[8] if de==d]
expected=[[1,d,i] for d,i in own if own[d,i][6]]
for D,V in [[4,4],[4,5]]:
 rr=[r for r in pairs if r[0]==D and r[2]==V];assert len(rr)==(630 if V==4 else 3240);assert [(r[1],r[3]) for r in rr]==[(i,j) for d,i in own if d==D for v,j in sm if v==V]
 for d,i,v,j,cl,yes in rr:
  actual=[a+b for a,b in zip(own[d,i][5],sm[v,j][5])];assert G.key(cl)==G.key(actual) and bool(yes)==(H.key(cl)==(0,0,0,0))
  if yes:expected.append([2,d,i,v,j])
assert expected==survivors and len(expected)==2944
ig,images,tests,remaining=data('filter_quotient14_isogeny_arithmetic.txt');iiv,irel,icd,icc,ihg,ihin=ig;assert iiv==iv and ihin==hinv;IG=Quotient(irel,[]);IH=Quotient(irel,ihg);assert IH.index==4 and all(IH.key(v)==(0,0,0,0) for v in icc);cache={}
for kk,e,w,cl in images:
 d=kk[0];row=[d,0,kk[1:]];E,T,Q=source(row);EP,(TP,U)=iso3(E,14*Q,[T,Q]);assert 2*TP==EP(0) and 14*U==EP(0) and 7*U!=EP(0) and 2*U!=EP(0) and TP!=7*U;de=degree(norm(EP,TP,U),d);assert e==de and w==d//e;cache[tuple(kk)]=cl
for k,desc,cl,yes in tests:
 assert desc==survivors[k-1];rows=[own[desc[1],desc[2]]]
 if desc[0]==2:rows.append(sm[desc[3],desc[4]])
 actual=[0]*len(iv)
 for row in rows:actual=[a+b for a,b in zip(actual,cache[tuple([row[0]]+row[2])])]
 assert IG.key(cl)==IG.key(actual) and bool(yes)==(IH.key(cl)==(0,0,0,0))
assert remaining==[r[1] for r in tests if r[3]] and len(remaining)==919
out=dict(complete=True,prime=5,target='X1(2,14)',target_genus=4,target_index=4,input_supports=11533,first_survivors=2944,isogeny_survivors=919,first_actual_maps=len(single)+len(small),independent_velu3_images=len(images))
if '--readonly' in sys.argv:assert json.loads((P/'filters_sage_certificate.json').read_text())==out
else:(P/'filters_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C242_LOWER_FILTERS_INDEPENDENT_SAGE_PASS',json.dumps(out),flush=True)
