from pathlib import Path
import ast,json,sys,time
from collections import Counter
from sage.all import GF,PolynomialRing,EllipticCurve,divisors,matrix,ZZ,vector
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
S=Path(__file__).resolve().parents[1]/'data';P=S/'cyclic21_filter';t0=time.monotonic()
def rd(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
class Quotient:
 def __init__(self,R,H):self.n=len(R[0])-1;self.M=hermite_normal_form(Matrix([r[:-1]for r in R+H]).T);self.index=abs(int(self.M.det()))
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0;w=list(v[:-1])
  for j in range(self.n-1,-1,-1):
   q=w[j]//int(self.M[j,j])
   for i in range(j+1):w[i]-=q*int(self.M[i,j])
  return tuple(w)
def add(*vs):return [sum(x)for x in zip(*vs)]
def scale(n,v):return [n*x for x in v]
def xy(E,Q):
 x,y=Q[0],Q[1];m=(3*x*x+2*E.a2()*x+E.a4())/(2*y);a2=3*x+E.a2()-m*m;a3=2*y;b=-a2**3/a3**2;c=1-2*m*a2/a3
 r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return[(s-r)/t,t/(s*s-s-r+1)]
def degree(vs,d):return min(int(e)for e in divisors(d)if all(v**(5**int(e))==v for v in vs))
r=rd(S/'filter_quotient14_linear_arithmetic.txt');own={(a[0],a[1]):a for a in r[2]};sm={(a[0],a[1]):a for a in r[3]};R=PolynomialRing(GF(5),'z');fields={d:GF(5**d,'z',modulus=R(pol))for d,pol in r[1]}
f=rd(P/'filter_arithmetic.txt');gd,images,tests,left=f;iv,rels,hg,cd,cc,comps=gd;assert iv==[2184,0]and Quotient(rels,hg).index==6
G=Quotient(rels,[]);H=Quotient(rels,hg);assert G.index==2184
qi,qr,qimages=rd(P/'quotient_coordinates.txt');QG=Quotient(qr,[]);assert QG.index==6;qm=matrix(ZZ,qimages).transpose()
def qvec(v):return list(qm*vector(ZZ,v))
assert all(QG.key(qvec(v))==(0,)for v in rels+hg)
Cp=[{H.key([0,0])}]+[set()for i in range(5)];CQ=[{QG.key([0,0])}]+[set()for i in range(5)]
for i,w in enumerate(cd):
 for d in range(w,6):
  for v in list(Cp[d-w]):Cp[d].add(H.key(add(list(v)+[0],cc[i])))
  for v in list(CQ[d-w]):CQ[d].add(QG.key(add(list(v)+[0],qvec(cc[i]))))
assert all(CQ[i]=={QG.key(v)for v in comps[i]}for i in range(6))
cache={};velu=0
for key,e,w,cl in images:
 mapno,d=key[:2];K=fields[d];z=K.gen();b=sum(K(v)*z**j for j,v in enumerate(key[2:2+d]));c=sum(K(v)*z**j for j,v in enumerate(key[2+d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P2=E(0,0);Q=E(b,b);assert Q.order()==42 and P2.order()==2 and P2!=21*Q;T3=2*Q
 if mapno==1:EG=E;U=T3
 else:
  TT=[P2,21*Q,P2+21*Q];assert len(set(TT))==3 and all(T.order()==2 for T in TT);T=TT[mapno-2];r=T[0];B=3*r*r+2*E.a2()*r+E.a4();x,y=T3[0],T3[1];v=x-r
  a4=E.a4()-5*B;a6=E.a6()-(4*E.a2()+7*r)*B;EG=EllipticCurve(K,[0,E.a2(),0,a4,a6]);U=EG(x+B/v,y*(1-B/(v*v)))
  iso=E.isogeny(T,check=True);assert iso.codomain()==EG and iso(T3)==U;velu+=1
 assert U.order()==21;de=degree(xy(EG,U),d);assert de==e and w==d//e;cache[tuple(key)]=cl
inputleft=rd(S/'genus5_quotient/filter_919/filter_arithmetic.txt')[-1];assert len(inputleft)==119 and len(tests)==119
for j,desc,cdeg,checks,passed in tests:
 assert desc==inputleft[j-1];rows=[own[desc[1],desc[2]]]
 if desc[0]==2:rows.append(sm[desc[3],desc[4]])
 assert cdeg==9-sum(a[0]for a in rows)
 for mapno,cl,yes in checks:
  actual=add(*(cache[tuple([mapno,a[0]]+a[2])]for a in rows));assert G.key(actual)==G.key(cl)and yes==(H.key(scale(-1,cl))in Cp[cdeg])
 assert [r[0]for r in checks]==list(range(1,len(checks)+1))and passed==all(r[2]for r in checks)and(passed or not checks[-1][2])
assert left==[r[1]for r in tests if r[-1]]and len(left)==44
allowed={tuple(d)for d in json.loads((S/'stabilizer_filter_sage.json').read_text())['remaining']};tail=[d for d in left if tuple(d)in allowed];assert len(tail)==18
out=dict(complete=True,prime=5,target='X1(21)',target_genus=5,finite_jacobian_order=2184,rational_image_order=364,index=6,input_supports=119,output_supports=44,after_optional_even_stabilizer=18,independent_normed_maps=len(images),independent_velu2_images=velu,cusp_complement_classes=[len(c)for c in Cp])
if '--readonly'in sys.argv:assert json.loads((P/'cyclic21_sage_certificate.json').read_text())==out
else:(P/'cyclic21_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n');(P/'after_cyclic21_and_stabilizer.txt').write_text(str(tail)+'\n')
print('C242_CYCLIC21_FOUR_MAP_INDEPENDENT_SAGE_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
