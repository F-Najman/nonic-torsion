"""Independent finite fields, full automorphism coverage, Velu2maps and HNF."""
from pathlib import Path
import ast,sys,json
from collections import Counter
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'common'))
from finite_arithmetic_char3 import Field,Curve,evaluate
P=Path(__file__).resolve().parent/'data'

def data(name):
 return ast.literal_eval((P/name).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def fkey(vs):
 out=[]
 for i in range(vs[0].K.n):out.append(tuple(z for v in vs for z in v.v));vs=[v**3 for v in vs]
 return min(out)
class Quotient:
 def __init__(self,R,hs):
  self.n=len(R[0])-1;self.H=hermite_normal_form(Matrix([r[:-1] for r in R+hs]).T);self.index=abs(int(self.H.det()))
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0;w=list(v[:-1])
  for j in range(self.n-1,-1,-1):
   q=w[j]//int(self.H[j,j])
   for i in range(j+1):w[i]-=q*int(self.H[i,j])
  return tuple(w)
 def add(self,u,v):return self.key([a+b for a,b in zip(u,v)]+[0])
def iso2(E,T,Q):
 K=E.K;r=T[0];A=E.a2+3*r;B=3*r*r+2*E.a2*r+E.a4
 xx=Q[0]-r;yy=Q[1];assert xx!=0;U=(xx+A+B/xx,yy*(1-B/(xx*xx)))
 EG=Curve(K,0,0);EG.a2=-2*A;EG.a4=A*A-4*B;EG.check(U)
 return EG,U

def main():
 inv,rels,hg,cd,cc,rs,remaining,coverage,fp=data('four_quotients_arithmetic.txt')
 assert inv==[4,4,24,936,0] and len(rs)==137 and len(coverage)==5
 QH=Quotient(rels,hg);assert QH.index==4
 allowed=[set() for _ in range(5)];allowed[0]={(0,)*QH.n}
 for w,v in zip(cd,cc):
  z=QH.key(v)
  for n in range(w,5):allowed[n]|={QH.add(x,z) for x in allowed[n-w]}
 assert [len(x) for x in allowed]==[1,1,1,3,3]
 assert [(r[0],r[1]) for r in remaining]==[(6,1),(7,10),(8,15),(8,24),(8,28),(8,29),(9,10),(9,19),(9,46),(9,66)]
 cov={r[0]:r for r in coverage};counts=[48,64,312,606,2112];orbits=[2,5,13,29,88]
 for d,n,o in zip(range(5,10),counts,orbits):
  _,pol,inputs,labels,os=cov[d];assert len(inputs)==n and len(set(map(tuple,inputs)))==n
  assert len(os)==o and set(labels)==set(range(1,o+1))
  assert {tuple(x) for orb in os for x in orb}==set(map(tuple,inputs))
 targets={};fields={};map_count=0;iso_count=0;laststab=None
 params=[[b,c,a] for b in range(2) for c in range(2) for a in [1,3,5,9,11,13]]
 for rec in rs:
  d,idx,pol,vals,maps,passed=rec
  if d not in fields:fields[d]=Field(pol)
  K=fields[d];assert K.pol==tuple(v%3 for v in pol) and K.n==d
  b,c=[evaluate(v,K.z) for v in vals];E=Curve(K,b,c);T=(K.zero,K.zero);U=(b,b);E.check(T);E.check(U)
  assert E.mul(2,T) is None and E.mul(28,U) is None and E.mul(14,U) is not None and E.mul(4,U) is not None
  R2=E.mul(14,U);assert R2!=T
  actual=[];stabs=[];orig=fkey([b,c])
  for bp,cp,a in params:
   T2=E.add(T,R2) if bp else T;UU=E.add(T,E.mul(a,U)) if cp else E.mul(a,U)
   xx=UU[0]-T2[0];yy=UU[1];iy=(yy*yy).inv();pair=[xx**3*iy,(c+3*T2[0])*xx*xx*iy]
   key=fkey(pair);actual.append(key)
   if key==orig:
    phases=[j for j in range(d) if pair==[b**(3**j),c**(3**j)]];assert len(phases)==1;stabs.append([bp,cp,a,phases[0]])
  _,pdef,inputs,labels,os=cov[d];assert pol==pdef
  assert set(actual)==set(map(tuple,os[idx-1]))
  assert set(actual)=={tuple(inputs[j]) for j in range(len(inputs)) if labels[j]==idx}
  if (d,idx)==(8,15):assert stabs==[[0,0,1,0],[1,0,13,4]];laststab=stabs
  variants=[(E,U),(E,E.add(T,U)),iso2(E,T,U),iso2(E,E.add(T,R2),U)];flags=[];assert len(maps)==4
  for variant,(EG,UU),row in zip(range(1,5),variants,maps):
   EG.check(UU);assert EG.mul(28,UU) is None and EG.mul(14,UU) is not None and EG.mul(4,UU) is not None
   pair=EG.tate_xy(UU);key=(d,)+fkey(pair);v,ti,e,mult,cl,flag=row;assert v==variant and d%e==0 and mult==d//e
   actuale=min(j for j in range(1,d+1) if d%j==0 and all(z**(3**j)==z for z in pair));assert e==actuale
   if ti in targets:assert targets[ti]==(key,e,cl)
   else:targets[ti]=(key,e,cl)
   accepted=QH.key([-z for z in cl]) in allowed[9-d];assert bool(flag)==accepted;flags.append(accepted)
   map_count+=1;iso_count+=variant>=3
  assert bool(passed)==all(flags)
  print('C228_D9_INDEPENDENT_SUPPORT_PASS',d,idx,flush=True)
 assert map_count==548 and iso_count==274 and laststab is not None
 traces=[];closed=[]
 for n in range(1,10):
  tr=-n*fp[50-n]-sum(fp[50-j]*traces[n-j-1] for j in range(1,n));traces.append(tr)
  num=3**n+1-tr-sum(d*closed[d-1] for d in range(1,n) if n%d==0);assert num%n==0;closed.append(num//n)
 cusps=Counter({1:12,2:6,3:4,6:2});nc=[v-cusps.get(d,0) for d,v in enumerate(closed,1)]
 assert nc==[0,0,0,0,48,64,312,606,2112]
 sigs,pairs,vs,tests,left=data('hecke_survivors_arithmetic.txt');assert len(sigs)==4 and sorted(map(len,pairs))==[2,2]
 expected=[]
 for i,j in pairs:
  for eps in [-1,1]:
   v=[0]*4;v[i-1]=-eps;v[j-1]=eps;expected.append(v)
 assert vs==expected;killed=set()
 for d,idx,q,corr,flags in tests:
  wanted=([[0]*4]+[[v[i]+q*w[i] for i in range(4)] for v in vs for w in vs]) if d==6 else [[0]*4]
  assert corr==wanted and len(flags)==len(wanted) and all(z in [0,1] for z in flags)
  if not any(flags):killed.add((d,idx))
 assert killed=={(r[0],r[1]) for r in remaining}-{(8,15)}
 assert [(r[0],r[1]) for r in left]==[(8,15)]
 out=dict(degree=9,complete=True,source_place_counts=nc,full_orbit_counts=orbits,supports=137,actual_maps=548,independent_two_isogenies=274,target_image_index=4,cusp_complement_class_counts=[1,1,1,3,3],distinct_target_keys=len(targets),Hecke_excluded_supports=9,corrected_stabilizer_excluded_support=[8,15],remaining=0,excluded_groups=[[2,28],[1,56]])
 if '--readonly' in sys.argv:assert json.loads((P/'replay_full.json').read_text())==out
 else:(P/'replay_full.json').write_text(json.dumps(out,indent=2)+'\n')
 print('C228_D9_FULL_FINITE_ARITHMETIC_REPLAY_PASS supports=137 orbit_keys=3142 actual_maps=548 two_isogenies=274')
 print('C228_D9_GROUP_AND_CUSP_COMPLEMENT_REPLAY_PASS target_index=4 remaining_after_Hecke_and_stabilizer=0')

if __name__=='__main__':main()
