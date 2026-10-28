from pathlib import Path
from math import gcd
from collections import Counter
import json
P=Path(__file__).resolve().parents[1]/'data'/'diagonal_quotient';CODE=Path(__file__).resolve().parent;n=36
labels={(a,b,c,d) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1 for c in range(2) for d in range(2) if (c*b-d*a)%2==1}
def orb(x):
 a,b,c,d=x;return frozenset(((s*(a+k*b))%n,s*b%n,s*(c+k*d)%2,s*d%2) for s in [-1,1] for k in range(n))
cs=sorted(set(map(orb,labels)),key=lambda x:sorted(x));ix={x:i for i,C in enumerate(cs) for x in C};assert len(cs)==64
h=[]
for C in cs:
 out={ix[(13*a%n,13*b%n,c,d)] for a,b,c,d in C};assert len(out)==1;h.append(out.pop())
assert all(h[h[h[i]]]==i for i in range(64))
left=set(range(64));qs=[]
while left:
 i=min(left);C={i,h[i],h[h[i]]};left-=C;qs.append(C)
assert len(qs)==32 and Counter(map(len,qs))=={1:16,3:16};qx={i:j for j,C in enumerate(qs) for i in C}
perms={}
for t in range(n):
 if gcd(t,n)!=1:continue
 pp=[]
 for C in qs:
  out={qx[ix[(t*a%n,b,t*c%2,d)]] for i in C for a,b,c,d in cs[i]};assert len(out)==1;pp.append(out.pop())
 assert len(set(pp))==32;perms[t]=pp
left=set(range(32));globalorbs=[]
while left:
 i=min(left);C={p[i] for p in perms.values()};left-=C;globalorbs.append(C)

from fractions import Fraction
from collections import defaultdict
def sig(x,i):
 a,b,c,d=x
 for w in range(1,37):
  y=((a+w*b)%36,b,(c+w*d)%2,d)
  if any(y==(s*a%36,s*b%36,s*c%2,s*d%2) for s in [-1,1]):break
 ram=3 if h[i]==i else 1;W=Fraction(w,ram);assert W.denominator==1
 def unitvals(N,r,qscale=Fraction(1)):
  r=min(r%N,N-r%N);lam=Fraction(r,N) if d==0 else Fraction(0);bv=Fraction(-2*r if d==0 else r,N);dv=-4*lam if d==0 else Fraction(1,2);dis=1-4*lam if d==0 else 0
  return [W*qscale*v for v in [bv,dv,dis,1-12*lam]]
 k=Fraction(9) if b%9==0 else Fraction(1) if b%3==0 else Fraction(1,9);r4=(b if b%3==0 and b%9 else 3*b)%4
 k3=Fraction(3) if b%3==0 else Fraction(1,3);r12p=(7*b//3 if b%3==0 else 7*b)%12
 return tuple(unitvals(12,b)+[W*(k-1)]+unitvals(4,r4,k)+unitvals(12,r12p,k3))
sigs=[]
for C in qs:
 vals={sig(x,i) for i in C for x in cs[i]};assert len(vals)==1,vals;sigs.append(vals.pop())

import ast,itertools
assert Counter(map(len,globalorbs))=={1:4,2:10,4:2}
assert len(set(sigs))==16
for C in globalorbs:assert len({sigs[i] for i in C})==1
raw=(P/'cusp_unit_signatures.txt').read_text().replace('[*','[').replace('*]',']');degrees,ns=ast.literal_eval(raw);ns=list(map(tuple,ns));assert len(ns)==20 and len(set(ns))==16
for v in sigs:assert all(x.denominator==1 for x in v)
sigs=[tuple(map(int,v)) for v in sigs]
left=set(range(32));cycles=[]
while left:
 i=min(left);C=set()
 while i not in C:C.add(i);left.remove(i);i=perms[5][i]
 cycles.append(C)
assert Counter(map(len,cycles))=={1:8,2:12};cix={i:j for j,C in enumerate(cycles) for i in C}
delta=[]
for C in qs:
 out={qx[ix[(7*a%36,7*b%36,c,d)]] for i in C for a,b,c,d in cs[i]};assert len(out)==1;delta.append(out.pop())
def closed_perm(pp):
 out=[]
 for C in cycles:
  target={pp[i] for i in C};j=cix[min(target)];assert target==cycles[j];out.append(j)
 return out
gamma=closed_perm(perms[7]);diam=closed_perm(delta)
tb=defaultdict(list);nb=defaultdict(list)
for i,C in enumerate(cycles):tb[sigs[min(C)],len(C)].append(i)
for i,(v,d) in enumerate(zip(ns,degrees)):nb[v,d].append(i)
assert set(tb)==set(nb) and all(len(tb[k])==len(nb[k]) for k in tb)
ks=sorted(tb);answers=set();assignments=0
for choice in itertools.product(*(list(itertools.permutations(nb[k])) for k in ks)):
 mp={i:j for k,perm in zip(ks,choice) for i,j in zip(tb[k],perm)};g=[None]*20;d=[None]*20
 for i,j in mp.items():g[j]=mp[gamma[i]];d[j]=mp[diam[i]]
 answers.add((tuple(g),tuple(d)));assignments+=1
assert len(answers)==1;g,d=answers.pop();rat=[i for i,C in enumerate(cycles) if len(C)==1 and all(pp[min(C)]==min(C) for pp in perms.values())];rats=sorted(mp[i]+1 for i in rat);assert len(rats)==4
assert all(g[g[i]]==i and d[d[i]]==i and g[d[i]]==d[g[i]] for i in range(20))
out=dict(complete=True,global_closed_cusps={'1':4,'2':10,'4':2},reduction5_closed_cusps={'1':8,'2':12},signature_functions=13,distinct_global_signatures=16,alignment_assignments_checked=assignments,alignment_independent=True,gamma7=[i+1 for i in g],diamond7=[i+1 for i in d],globally_rational_cusp_indices=rats)
import sys
permtext='CuspGamma7:='+str(out['gamma7'])+';\nCuspDiamond7:='+str(out['diamond7'])+';\nCuspGlobalRational:='+str(rats)+';\n'
if '--readonly' in sys.argv:
 assert json.loads((P/'tate_cusp_alignment.json').read_text())==out and (CODE/'cusp_permutations.m').read_text()==permtext
else:
 (P/'tate_cusp_alignment.json').write_text(json.dumps(out,indent=2)+'\n');(CODE/'cusp_permutations.m').write_text(permtext)
print('C236_DELTA13_TATE_CUSP_UNIT_ALIGNMENT_PASS',json.dumps(out),flush=True)
