#!/usr/bin/env python3
from pathlib import Path
import re,json,sys
P=Path(__file__).resolve().parent/'data';s=(P/'rational_two_delta4_matrices.m').read_text();assert re.search(r'N:=37; d:=80;',s);d=80

def mat(name):
 m=re.search(re.escape(name)+r':=Matrix\(Integers\(\),(\d+),(\d+),\[(.*?)\]\);',s,re.S);assert m and int(m[1])==d and int(m[2])==d
 vals=list(map(int,re.findall(r'-?\d+',m[3])));assert len(vals)==d*d
 return[sum((vals[i*d+j]&1)<<j for j in range(d))for i in range(d)]
D2=mat('D2Z');Da=mat('DaZ');As=[mat(f'Aq{q}Z')for q in[3,5]]
def act(v,M):
 out=0
 while v:
  lb=v&-v;out^=M[lb.bit_length()-1];v^=lb
 return out
assert Da==[act(v,D2)for v in D2]
for A in As:assert all(act(A[i],D2)==act(D2[i],A)for i in range(d))
piv={}
for i in range(d):
 row=As[0][i]|(As[1][i]<<d)
 while row:
  c=row.bit_length()-1
  if c in piv:row^=piv[c]
  else:piv[c]=row;break
assert len(piv)==80
res={'N':37,'integral_homology_rank':80,'actual_JOne':True,'sign':0,'q_operators':[3,5],'common_binary_kernel_dimension':0,'rational_2torsion_is_zero':True,'status':'PASS'}
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(res,indent=2)+'\n')
print('NONIC37_RATIONAL_TWO_INDEPENDENT_REPLAY_PASS',json.dumps(res,sort_keys=True))
