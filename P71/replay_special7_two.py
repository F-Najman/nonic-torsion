#!/usr/bin/env python3
"""Independent binary replay on the exported actual integral JOne homology."""
from pathlib import Path
import re,sys,json
D=Path(__file__).resolve().parent/'data';N=71
s=(D/'special7_two_matrices.m').read_text();d=int(re.search(r'\bd:=(\d+);',s).group(1))
assert d==352
def mat(name):
 m=re.search(re.escape(name)+r':=Matrix\(Integers\(\),(\d+),(\d+),\[(.*?)\]\);',s,re.S);assert m,name
 assert int(m[1])==d and int(m[2])==d
 vals=list(map(int,re.findall(r'-?\d+',m[3])));assert len(vals)==d*d
 return[sum((vals[i*d+j]&1)<<j for j in range(d))for i in range(d)]
D2=mat('D2Z');Da=mat('DaZ')
qs=[int(q) for q in re.findall(r'Aq(\d+)Z:=',s)];assert qs==[3];As=[mat(f'Aq{q}Z')for q in qs]
def act(v,M):
 out=0
 while v:
  lb=v&-v;out^=M[lb.bit_length()-1];v^=lb
 return out
for A in As:assert all(act(A[i],D2)==act(D2[i],A)for i in range(d))
power=[1<<i for i in range(d)]
for _ in range(5):power=[act(v,D2)for v in power]
assert power==Da
r=7;power=[1<<i for i in range(d)]
for _ in range(r):power=[act(v,Da)for v in power]
assert power==[1<<i for i in range(d)]
piv={};ker=[]
for i in range(d):
 row=As[0][i];comb=1<<i
 while row:
  c=row.bit_length()-1
  if c in piv:old,rel=piv[c];row^=old;comb^=rel
  else:piv[c]=(row,comb);break
 if row==0:ker.append(comb)
assert len(ker)==0
assert all(act(v,[D2[i]^(1<<i)for i in range(d)])==0 for v in ker)
res={'N':N,'integral_homology_rank':d,'actual_JOne':True,'sign':0,'q_operators':qs,'common_binary_kernel_dimension':0,'diamond32_difference_rank_on_kernel':0,'diamond_order':r,'status':'PASS'}
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(res,indent=2)+'\n')
print('NONIC_RATIONAL_TWO_INDEPENDENT_REPLAY_PASS',json.dumps(res,sort_keys=True))
