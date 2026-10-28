"""Independent Tate-cusp labels, bvaluations, globalQalignment and hcycles.

For Q=zeta26^a q^(b/26), P2=(-1)^c q^(d/2), put r=min(b,26-b).
The normalized modularunit b has valuation(w/26)*(-2r if d=0 else r).
This follows from xTate(Q)-xTate(P2) and yTate(Q)+xTate(Q)/2; all leading
constants are units awayfrom2,13. Widthw is the unipotent orbitlength.
"""
from pathlib import Path
from math import gcd
from fractions import Fraction
from collections import Counter
import ast,json,sys
P=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/'data'/'occurrence'

def data():
 s=(P/'finite_geometric_arithmetic.txt').read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False');return ast.literal_eval(s)

n=26
labels={(a,b,c,d) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1 for c in range(2) for d in range(2) if (c*b-d*a)%2==1}
def orbit(x):
 a,b,c,d=x
 return frozenset(((s*(a+k*b))%n,s*b%n,s*(c+k*d)%2,s*d%2) for s in (-1,1) for k in range(n))
cs=sorted(set(map(orbit,labels)),key=lambda x:sorted(x));index={x:i for i,C in enumerate(cs) for x in C}
assert len(cs)==36

def val_b(x):
 a,b,c,d=x;r=min(b,n-b)
 for w in range(1,n+1):
  z=((a+w*b)%n,b,(c+w*d)%2,d)
  if any(z==(s*a%n,s*b%n,s*c%2,s*d%2) for s in (-1,1)):break
 v=Fraction(w,n)*(-2*r if d==0 else r);assert v.denominator==1
 return int(v),w
vv=[];width=[]
for C in cs:
 out={val_b(x) for x in C};assert len(out)==1;v,w=out.pop();vv.append(v);width.append(w)
assert sum(width)==504 and sum(max(v,0) for v in vv)==84 and sum(min(v,0) for v in vv)==-84

def gal(t):
 out=[]
 for C in cs:
  ans={index[(t*a%n,b,t*c%2,d)] for a,b,c,d in C};assert len(ans)==1;out.append(ans.pop())
 return out
perms={t:gal(t) for t in range(n) if gcd(t,n)==1}
rat=[i for i in range(36) if all(pp[i]==i for pp in perms.values())];assert len(rat)==18
assert len({vv[i] for i in rat})==18
left=set(range(36));gs=[]
while left:
 i=min(left);C={pp[i] for pp in perms.values()};left-=C;gs.append(C)
assert Counter(map(len,gs))=={1:18,6:3}
for p,sig in [(3,{1:18,3:6}),(7,{1:18,6:3})]:
 left=set(range(36));cycles=[];pp=perms[p]
 while left:
  i=min(left);C=[]
  while i not in C:C.append(i);left.remove(i);i=pp[i]
  cycles.append(C)
 assert Counter(map(len,cycles))==sig
 assert {C[0] for C in cycles if len(C)==1}==set(rat)
h=[]
for C in cs:
 ans={index[((13*c+3*a)%26,(13*d+3*b)%26,a%2,b%2)] for a,b,c,d in C};assert len(ans)==1;h.append(ans.pop())
assert all(h[h[h[i]]]==i for i in range(36)) and all(h[i]!=i for i in range(36))
byval={vv[i]:i for i in rat};H={v:vv[h[i]] for v,i in byval.items()}
ls,av,am,bv,cv,cm,records=data();assert ls==sorted(byval)
assert records[0][0]==3 and records[0][-1]==[H[v] for v in ls]
A=dict(zip(av,am));B={v:1 for v in bv};C=dict(zip(cv,cm))
assert sum(A.values())==12 and sum(B.values())==3 and sum(C.values())==12 and A!=C
assert all(A.get(H[v],0)==A.get(v,0) and B.get(H[v],0)==B.get(v,0) and C.get(H[v],0)==C.get(v,0) for v in byval)
assert len({bv[0],H[bv[0]],H[H[bv[0]]]})==3 and {bv[0],H[bv[0]],H[H[bv[0]]]}==set(bv)
assert not (set(A)&set(B)) and not (set(C)&set(B))
print('C226_TATE_GLOBAL_CUSP_ALIGNMENT_PASS rational18 nonrational3x6 goodprimes3,7 bdegree84')
print('C226_TATE_SYMMETRIC_PENCIL_PASS h_order3=true A_stable=true B_single3cycle=true C_stable=true')
print(json.dumps(dict(b_labels=ls,h_action=H,cusp_width_sum=sum(width),geometric_cusps=36),indent=2))
