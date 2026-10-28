"""Independent Sage replay of tail_all267.m (rerun made for the repository, 27 October 2026).

Adapted from replay_cyclic_sage.py, which replays the two runs of the record (the interrupted
prefix of filter_cyclic36_linear.m and tail.m); here the single stage is the run of tail_all267.m
on all 267 genus-9 survivors of the diagonal quotient filter.  Checks, without Magma: every image
(Tate normal form, 2-isogenies by Velu's formulas, field of definition), the class sums, the
membership of minus each class in the cusp complements modulo the rational cuspidal subgroup
(Hermite normal form), and that exactly one of the 267 survivors, [1, 8, 1601], passes all four
maps.  Run from the repository root: python3 C2x36/cyclic_tail_after_genus9/replay_cyclic_all267_sage.py
"""
from pathlib import Path
import ast,json,sys
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'C2x36'/'data'

def data(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
class Quotient:
 def __init__(self,R,hs):
  self.n=len(R[0])-1;self.H=hermite_normal_form(Matrix([r[:-1] for r in R+hs]).T);self.index=abs(int(self.H.det()))
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0;w=list(v[:-1])
  for j in range(self.n-1,-1,-1):
   q=w[j]//int(self.H[j,j])
   for i in range(j+1):w[i]-=q*int(self.H[i,j])
  return tuple(w)
def tate(E,Q):
 x,y=Q[0],Q[1];a2,a4=E.a2(),E.a4();aa3=2*y;aa2=3*x+a2;aa4=3*x*x+2*x*a2+a4;aa1=2*aa4/aa3;aa2-=(aa4/aa3)**2
 b=-aa2**3/aa3**2;c=-(aa1*aa2-aa3)/aa3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1
 return [(s-r)/t,t/(s*s-s-r+1)]
def iso2(E,T,Q):
 K=E.base_field();r=T[0];A=E.a2()+3*r;B=3*r*r+2*E.a2()*r+E.a4();x=Q[0]-r;u=(x+A+B/x,Q[1]*(1-B/x**2));EP=EllipticCurve(K,[0,-2*A,0,A*A-4*B,0]);return EP,EP(u)
first=data(D/'filter_quotient18_linear_arithmetic.txt');own={(r[0],r[1]):r for r in first[2]};small={(r[0],r[1]):r for r in first[3]};R=PolynomialRing(GF(5),'z');fields={d:GF(5**d,'z',modulus=R(f)) for d,f in first[1]}
projected=data(D/'diagonal_quotient'/'quotient_filter_arithmetic.txt')[3];assert len(projected)==267
stages=[];alltested={};surv=[]
for label,filename,expected in [('all267',D/'cyclic_tail_after_genus9'/'tail_all267_arithmetic.txt',projected)]:
 group,images,tests,remaining=data(filename);iv,rels,hg,cd,cc,nativecomps=group;assert iv==[4,4,12,12,36,7812,7812,0]
 G=Quotient(rels,[]);H=Quotient(rels,hg);assert H.index==214272
 zero=(0,)*7;comp=[{zero},set(),set(),set()]
 for degree,vec in zip(cd,cc):
  if degree>3:continue
  for d in range(degree,4):
   comp[d].update(H.key([a+b for a,b in zip(list(z)+[0],vec)]) for z in list(comp[d-degree]))
 assert [len(x) for x in comp]==[1,7,27,77]
 cached={};isogs=0
 for key,e,w,cl in images:
  v,d=key[:2];ky=key[2:];K=fields[d];z=K.gen();b=sum(K(k)*z**j for j,k in enumerate(ky[:d]));c=sum(K(k)*z**j for j,k in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b)
  if v==1:EP,U=E,Q
  elif v==2:EP,U=E,T+Q
  elif v in [3,4]:EP,U=iso2(E,T if v==3 else T+18*Q,Q);isogs+=1
  else:raise AssertionError(v)
  assert 36*U==EP(0) and 18*U!=EP(0) and 12*U!=EP(0);xy=tate(EP,U);de=min(int(j) for j in divisors(d) if all(t**(5**int(j))==t for t in xy));assert de==e and w==d//e
  cached[tuple(key)]=cl
 for k,desc,cd,checks,yes in tests:
  assert desc==expected[k-1];rows=[own[desc[1],desc[2]]]
  if desc[0]>=2:rows.append(small[desc[3],desc[4]])
  if desc[0]==3:rows.append(small[desc[5],desc[6]])
  assert cd==9-sum(r[0] for r in rows);assert [q[0] for q in checks]==list(range(1,len(checks)+1))
  for v,cl,flag in checks:
   summed=[0]*len(iv)
   for row in rows:summed=[x+y for x,y in zip(summed,cached[tuple([v,row[0]]+row[2])])]
   assert G.key(cl)==G.key(summed);assert bool(flag)==(H.key([-x for x in cl]) in comp[cd])
  assert bool(yes)==all(q[2] for q in checks) and (len(checks)==4 or not yes);alltested[tuple(desc)]=bool(yes)
 assert remaining==[r[1] for r in tests if r[4]];surv+=remaining
 stages.append(dict(stage=label,tested=len(tests),expected=len(expected),complete=len(tests)==len(expected),remaining=remaining,actual_images=len(images),independent_velu2_images=isogs,index=H.index,cusp_complement_counts=[len(x) for x in comp]))
 print('C236_CYCLIC_SAGE_STAGE_PASS',stages[-1],flush=True)
assert all(s['complete'] for s in stages)
assert len(alltested)==267 and set(map(tuple,projected))==set(alltested)
assert [d for d in projected if alltested[tuple(d)]]==[[1,8,1601]]
out=dict(complete=True,stages=stages,genus9_survivors=267,cyclic_survivors=[[1,8,1601]],covered_genus9_survivors=267)
cert=D/'cyclic_tail_after_genus9'/'cyclic_all267_sage_certificate.json'
if '--readonly' in sys.argv:assert json.loads(cert.read_text())==out
else:cert.write_text(json.dumps(out,indent=2)+'\n')
print('C236_CYCLIC_SAGE_ALL267_REPLAY_PASS',json.dumps(out),flush=True)
