#!/usr/bin/env python3
"""Read-only exact q certificate for the degree9 canonical j-map."""
from pathlib import Path
import ast,re,json,hashlib
H=Path(__file__).resolve().parents[1]/'data'/'quotient19';cs=(H/'canonical_data.m').read_text();js=(H/'jmap_q_data.m').read_text();N=600;mask=(1<<N)-1

def lit(s,n):return ast.literal_eval(re.search(re.escape(n)+r':=(.*?);',s,re.S)[1].replace('<','(').replace('>',')'))
def matrix(s,n):
 m=re.search(re.escape(n)+r':=Matrix\((?:GF\(2\)|Integers\(\)),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,n
 r,c=map(int,m.group(1,2));a=ast.literal_eval(m[3]);assert len(a)==r*c
 return r,c,[sum((a[i*c+j]&1)<<j for j in range(c))for i in range(r)]
def rank(rows):
 p={}
 for x in rows:
  while x:
   k=x.bit_length()-1
   if k in p:x^=p[k]
   else:p[k]=x;break
 return len(p)
def combine(v,rows):
 z=0
 for a,b in zip(v,rows):
  if a:z^=b
 return z

def mul(a,b):
 if a.bit_count()<b.bit_count():a,b=b,a
 z=0
 while b:
  k=b&-b;z^=a<<(k.bit_length()-1);b^=k
 return z&mask
old=matrix(cs,'qcoeffs');q=matrix(js,'qcoeffs600');md=matrix(js,'Diamond5Matrix');rel=matrix(cs,'QuadRelations');assert old[:2]==(19,130)and q[:2]==(19,600)and md[:2]==(19,19)and rel[:2]==(136,190)
assert [x&((1<<130)-1)for x in q[2]]==old[2]and rank([x&((1<<50)-1)for x in old[2]])==19
I=[1<<i for i in range(19)];p=I
for n in range(1,5):
 p=[combine([(x>>j)&1 for j in range(19)],md[2])for x in p]
 if n<4:assert p!=I
assert p==I
pairs=lit(cs,'pairs');assert len(pairs)==190 and set(pairs)=={(i,j)for i in range(1,20)for j in range(i,20)}
products=[mul(q[2][i-1],q[2][j-1])for i,j in pairs];assert rank(products)==54
for row in rel[2]:assert combine([(row>>j)&1 for j in range(190)],products)==0
powers=[[1]for _ in range(19)]
for i in range(19):
 for k in range(1,10):powers[i].append(mul(powers[i][-1],q[2][i]))
def evalmon(e):
 assert len(e)==19;x=1
 for i,k in enumerate(e):
  if k:x=mul(x,powers[i][k])
 return x
E3=list(map(evalmon,lit(js,'Std3Exponents')));E9=list(map(evalmon,lit(js,'Std9Exponents')));assert len(E3)==90 and len(E9)==306 and rank(E3)==90 and rank(E9)==306
Delta=sum(1<<(i*i)for i in range(1,25,2)if i*i<N);DE3=[mul(Delta,x)for x in E3]
assert len(E3)+2*len(E9)-rank([a|(b<<600)for a,b in zip(E3,DE3)]+E9+[x<<600 for x in E9])==10
G=combine(lit(js,'Gvec'),E3);A=combine(lit(js,'NumeratorVec'),E9);B=combine(lit(js,'DenominatorVec'),E9);assert G and A==G and B==mul(Delta,G)and B
print('NONIC73_Q_IDENTITIES_READONLY_PASS canonical_degree9 auxiliary_degree3 simultaneous_dimension10 q_precision600 Sturm444 quadrics136 diamond5_order4')
