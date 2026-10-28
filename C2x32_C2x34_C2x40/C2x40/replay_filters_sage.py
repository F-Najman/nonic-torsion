from pathlib import Path
import ast,json,sys
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parents[1]/'data'/'C2x40'
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
def iso2(E,T,U):
 K=E.base_field();r=T[0];A=E.a2()+3*r;B=3*r*r+2*E.a2()*r+E.a4();EP=EllipticCurve(K,[0,-2*A,0,A*A-4*B,0])
 def point(Q):x=Q[0]-r;return EP(x+A+B/x,Q[1]*(1-B/x**2))
 return EP,[point(Q) for Q in U]
f=data('filter_quotient20_linear_arithmetic.txt');g,fs,single,small,pairs,triples,survivors,coverage=f;iv,rels,cd,cc,hg,hinv=g;assert iv==[8,120,120,0] and hinv==[4,60,120] and triples==[];H=Quotient(rels,hg);G=Quotient(rels,[]);assert H.index==4 and all(H.key(v)==(0,0,0) for v in cc)
R=PolynomialRing(GF(3),'z');fields={d:GF(3**d,'z',modulus=R(pol)) for d,pol in fs};own={};sm={}
def source(row):
 d,i,ky=row[:3];K=fields[d];z=K.gen();b=sum(K(v)*z**j for j,v in enumerate(ky[:d]));c=sum(K(v)*z**j for j,v in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);return E,E(0,0),E(b,b)
def degree(pair,d):return min(int(e) for e in divisors(d) if all(v**(3**int(e))==v for v in pair))
for is_small,rows in [(False,single),(True,small)]:
 for row in rows:
  d,i,ky,e,w,cl=row[:6];E,T,Q=source(row);de=degree(norm(E,T,2*Q),d);assert e==de and w==d//e;(sm if is_small else own)[d,i]=row
  if not is_small:assert bool(row[6])==(H.key(cl)==(0,0,0))
expected=[[1,d,i] for d,i in own if own[d,i][6]]
for D,V in [[4,4],[4,5]]:
 rr=[r for r in pairs if r[0]==D and r[2]==V];assert len(rr)==(32 if V==4 else 96);assert [(r[1],r[3]) for r in rr]==[(i,j) for d,i in own if d==D for v,j in sm if v==V]
 for d,i,v,j,cl,yes in rr:
  actual=[a+b for a,b in zip(own[d,i][5],sm[v,j][5])];assert G.key(cl)==G.key(actual) and bool(yes)==(H.key(cl)==(0,0,0))
  if yes:expected.append([2,d,i,v,j])
assert expected==survivors and len(expected)==64
ig,images,tests,remaining=data('filter_quotient20_isogeny_arithmetic.txt');iiv,irel,icd,icc,ihg,ihin=ig;assert iiv==iv and ihin==hinv;IG=Quotient(irel,[]);IH=Quotient(irel,ihg);assert IH.index==4 and all(IH.key(v)==(0,0,0) for v in icc);cache={}
for kk,e,w,cl in images:
 d=kk[0];row=[d,0,kk[1:]];E,T,Q=source(row);EP,(TP,U)=iso2(E,20*Q,[T,Q]);assert 2*TP==EP(0) and 20*U==EP(0) and 10*U!=EP(0) and 4*U!=EP(0) and TP!=10*U;de=degree(norm(EP,TP,U),d);assert e==de and w==d//e;cache[tuple(kk)]=cl
for k,desc,cl,yes in tests:
 assert desc==survivors[k-1];rows=[own[desc[1],desc[2]]]
 if desc[0]==2:rows.append(sm[desc[3],desc[4]])
 actual=[0]*len(iv)
 for row in rows:actual=[a+b for a,b in zip(actual,cache[tuple([row[0]]+row[2])])]
 assert IG.key(cl)==IG.key(actual) and bool(yes)==(IH.key(cl)==(0,0,0))
assert remaining==[r[1] for r in tests if r[3]] and len(remaining)==27
out=dict(complete=True,prime=3,target='X1(2,20)',target_genus=9,target_index=4,input_supports=234,first_survivors=64,isogeny_survivors=27,first_actual_maps=len(single)+len(small),independent_velu2_images=len(images))
if '--readonly' in sys.argv:assert json.loads((P/'filters_sage_certificate.json').read_text())==out
else:(P/'filters_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C240_LOWER_FILTERS_INDEPENDENT_SAGE_PASS',json.dumps(out),flush=True)
