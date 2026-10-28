#!/usr/bin/env python3
"""Independent binary replay of all12 RR kernels and20 cusp restrictions."""
import itertools,json
from pathlib import Path
H=Path(__file__).resolve().parent;D=json.loads((H.parent/'data'/'C95'/'fi3_data.json').read_text())
def rank(A):
 A=[r[:] for r in A];k=0
 for j in range(len(A[0]) if A else 0):
  i=next((i for i in range(k,len(A)) if A[i][j]),None)
  if i is None:continue
  A[i],A[k]=A[k],A[i]
  for i in range(k+1,len(A)):
   if A[i][j]:A[i]=[x^y for x,y in zip(A[i],A[k])]
  k+=1
 return k
def matrix(v,n,m):assert len(v)==n*m;return [v[i*m:(i+1)*m] for i in range(n)]
def mul(A,B):return [[sum(x*y for x,y in zip(r,c))%2 for c in zip(*B)] for r in A]
C=matrix(D['Cartier'],9,9);U=matrix(D['U_basis'],D['U_dimension'],9);CI=[[C[i][j]^(i==j) for j in range(9)] for i in range(9)]
assert rank(C)==7 and rank(CI)==rank(U)==8 and rank(CI+U)==8
assert sorted(D['cusp_widths'])==[1,5,19,95]
blocks={}
for e in D['kernels']:
 i,m=e['cusp'],e['multiplicity'];W=matrix(e['kernel'],e['kernel_dimension'],9);A=matrix(e['restriction'],9,m)
 assert e['kernel_dimension']==9-m and rank(W)==9-m and rank(A)==m and not any(any(r) for r in mul(W,A))
 blocks[i,m]=mul(U,A)
assert len(blocks)==12
expected={v for v in itertools.product(range(4),repeat=4) if sum(v)==3};assert len(expected)==20
assert {tuple(r['multiplicities']) for r in D['configurations']}==expected
for r in D['configurations']:
 cols=[blocks[i,m] for i,m in enumerate(r['multiplicities'],1) if m]
 A=[sum((B[j] for B in cols),[]) for j in range(8)]
 assert rank(A)==r['rank']==3
print('NONIC_C95_BINARY_FI_PASS Cartier_rank=7 U_dimension=8 RR_kernels=12 configurations=20 all_restriction_ranks=3 complete=true')
