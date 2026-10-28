#!/usr/bin/env python3
"""Independent exact basis alignment, old/new RR transport and3098 FI ranks."""
import ast,itertools,json,re
from pathlib import Path
H=Path(__file__).resolve().parent;D=json.loads((H/'data'/'fi9_data.json').read_text());src=(H/'data'/'retained_degree_independent_kernels.m').read_text();n=37

def matrix(v,r,c):assert len(v)==r*c;return [v[i*c:(i+1)*c] for i in range(r)]
def parse(name):
 m=re.search(re.escape(name)+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',src,re.S);assert m,name
 return matrix(ast.literal_eval(m[3]),int(m[1]),int(m[2]))
def mul(A,B):return [[sum(x*y for x,y in zip(r,c))%2 for c in zip(*B)] for r in A]
def add(A,B):return [[x^y for x,y in zip(r,s)] for r,s in zip(A,B)]
def rank(A):
 piv={}
 for r in A:
  v=sum(x<<j for j,x in enumerate(r))
  while v:
   j=v.bit_length()-1
   if j in piv:v^=piv[j]
   else:piv[j]=v;break
 return len(piv)
def brank(rows):
 piv={}
 for v in rows:
  while v:
   j=v.bit_length()-1
   if j in piv:v^=piv[j]
   else:piv[j]=v;break
 return len(piv)
def inv(A):
 B=[r[:]+[int(i==j) for j in range(n)] for i,r in enumerate(A)]
 for j in range(n):
  k=next(i for i in range(j,n) if B[i][j]);B[k],B[j]=B[j],B[k]
  for i in range(n):
   if i!=j and B[i][j]:B[i]=[x^y for x,y in zip(B[i],B[j])]
 return [r[n:] for r in B]
def power(A,k):
 B=[[int(i==j) for j in range(n)] for i in range(n)]
 for _ in range(k):B=mul(B,A)
 return B
I=[[int(i==j) for j in range(n)] for i in range(n)]
C=matrix(D['Cartier'],n,n);U=matrix(D['U_basis'],32,n);T=matrix(D['T_old_to_new'],n,n);M2=matrix(D['deck2'],n,n);M3=matrix(D['deck3'],n,n)
assert rank(U)==32 and rank(mul(C,add(C,I))+U)==32 and rank(T)==37
assert power(M2,2)==I and power(M3,3)==I and mul(M2,M3)==mul(M3,M2)
Cold=parse('CartierMatrix');O2=parse('MatrixC2');O3=parse('MatrixC3');W8=parse('sourcekernels[8]');An8=matrix(D['anchor_restriction'],n,8);Wn8=matrix(D['anchor_kernel'],29,n)
assert rank(W8)==rank(Wn8)==29 and rank(An8)==8 and not any(any(r) for r in mul(mul(W8,T),An8))
eq=[]
for A,B in [(Cold,C),(O2,M2),(O3,M3)]:
 assert mul(A,T)==mul(T,B)
 for i in range(n):
  for j in range(n):
   mask=0
   for k in range(n):
    if A[i][k]:mask^=1<<(k*n+j)
    if B[k][j]:mask^=1<<(i*n+k)
   eq.append(mask)
for r in W8:
 for c in range(8):
  cmask=sum(1<<j for j in range(n) if An8[j][c]);mask=0
  for i in range(n):
   if r[i]:mask^=cmask<<(i*n)
  eq.append(mask)
assert brank(eq)==n*n-1
Ti=inv(T);blocks={};As={}
for e in D['blocks']:
 i,m=e['cusp'],e['multiplicity'];length=m*D['degrees'][i-1];A=matrix(e['restriction'],n,length);W=matrix(e['kernel'],e['kernel_dimension'],n)
 assert rank(A)==length and rank(W)==e['kernel_dimension']==n-length
 assert not any(any(r) for r in mul(W,A))
 As[i,m]=A;blocks[i,m]=mul(U,A)
 if i<=6 and m<=8:
  oldA=[r[:m] for r in parse(f'basisjets[{i}]')];assert A==mul(Ti,oldA)
assert len(As)==60
p2=[2,1,4,3,6,5];p3=[3,4,5,6,1,2]
for i in range(1,7):
 matches=[]
 for a,b in itertools.product(range(2),range(3)):
  j=1
  for _ in range(a):j=p2[j-1]
  for _ in range(b):j=p3[j-1]
  if j==i:matches.append((a,b))
 assert len(matches)==1;a,b=matches[0];assert As[i,9]==mul(mul(power(M2,a),power(M3,b)),As[1,9])
for m in [1,2,3]:assert As[8,m]==mul(M2,As[7,m])
def weighted(ds,d):
 if not ds:
  if d==0:yield ()
  return
 for k in range(d//ds[0]+1):
  for rest in weighted(ds[1:],d-k*ds[0]):yield (k,)+rest
vectors=list(weighted([1]*6+[3,3],9));assert len(vectors)==D['tested']==3098 and D['exceptions']==[]
for v in vectors:
 B=[[] for _ in range(32)]
 for i,m in enumerate(v,1):
  if m:
   for j in range(32):B[j]+=blocks[i,m][j]
 assert len(B[0])==9 and rank(B)==9,v
print('NONIC77_INDEPENDENT_FI9_PASS unique_basis_intertwiner_dimension=1 RR_blocks=60 cubic_cusps=2 tested=3098 all_ranks=9 complete=true')
