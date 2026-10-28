#!/usr/bin/env python3
"""Independent binary replay on the exported actual integral JOne homology."""
from pathlib import Path
import re,sys,json
D=Path(__file__).resolve().parents[1];N=int(sys.argv[1]);assert N in[43,73]
s=(D/{43:'P43_P53_P61/data/P43/N43_matrices.m',73:'P73/data/deck/N73_matrices.m'}[N]).read_text();d=int(re.search(r'\bd:=(\d+);',s).group(1))
assert d=={43:114,73:374}[N]
def mat(name):
 m=re.search(re.escape(name)+r':=Matrix\(Integers\(\),(\d+),(\d+),\[(.*?)\]\);',s,re.S);assert m,name
 assert int(m[1])==d and int(m[2])==d
 vals=list(map(int,re.findall(r'-?\d+',m[3])));assert len(vals)==d*d
 return[sum((vals[i*d+j]&1)<<j for j in range(d))for i in range(d)]
D2=mat('D2Z');Da=mat('DaZ');assert D2==Da
As=[mat(f'Aq{q}Z')for q in[3,5]]
def act(v,M):
 out=0
 while v:
  lb=v&-v;out^=M[lb.bit_length()-1];v^=lb
 return out
for A in As:assert all(act(A[i],D2)==act(D2[i],A)for i in range(d))
r={43:7,73:9}[N];power=[1<<i for i in range(d)]
for _ in range(r):power=[act(v,D2)for v in power]
assert power==[1<<i for i in range(d)]
piv={};ker=[]
for i in range(d):
 row=As[0][i]|(As[1][i]<<d);comb=1<<i
 while row:
  c=row.bit_length()-1
  if c in piv:old,rel=piv[c];row^=old;comb^=rel
  else:piv[c]=(row,comb);break
 if row==0:ker.append(comb)
assert len(ker)==6
assert all(act(v,[D2[i]^(1<<i)for i in range(d)])==0 for v in ker)
res={'N':N,'integral_homology_rank':d,'actual_JOne':True,'sign':0,'q_operators':[3,5],'common_binary_kernel_dimension':6,'diamond2_difference_rank_on_kernel':0,'diamond_order':r,'status':'PASS'}
if len(sys.argv)>2:Path(sys.argv[2]).write_text(json.dumps(res,indent=2)+'\n')
print('NONIC_RATIONAL_TWO_INDEPENDENT_REPLAY_PASS',json.dumps(res,sort_keys=True))
