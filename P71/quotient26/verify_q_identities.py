#!/usr/bin/env python3
"""Read-only exact q certificate for the degree9 canonical j-map."""
from pathlib import Path
import ast,re,json,hashlib
H=Path(__file__).resolve().parents[1]/'data'/'quotient26';cs=(H/'canonical_data.m').read_text();js=(H/'jmap_q_data.m').read_text();N=600;mask=(1<<N)-1

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
old=matrix(cs,'qcoeffs');q=matrix(js,'qcoeffs600');md=matrix(js,'Diamond7Matrix');rel=matrix(cs,'QuadRelations');assert old[:2]==(26,130)and q[:2]==(26,600)and md[:2]==(26,26)and rel[:2]==(276,351)
assert [x&((1<<130)-1)for x in q[2]]==old[2]and rank([x&((1<<61)-1)for x in old[2]])==26
I=[1<<i for i in range(26)];p=I
for n in range(1,6):
 p=[combine([(x>>j)&1 for j in range(26)],md[2])for x in p]
 if n<5:assert p!=I
assert p==I
pairs=lit(cs,'pairs');assert len(pairs)==351 and set(pairs)=={(i,j)for i in range(1,27)for j in range(i,27)}
products=[mul(q[2][i-1],q[2][j-1])for i,j in pairs];assert rank(products)==75
for row in rel[2]:assert combine([(row>>j)&1 for j in range(351)],products)==0
powers=[[1]for _ in range(26)]
for i in range(26):
 for k in range(1,9):powers[i].append(mul(powers[i][-1],q[2][i]))
def evalmon(e):
 assert len(e)==26;x=1
 for i,k in enumerate(e):
  if k:x=mul(x,powers[i][k])
 return x
E2=list(map(evalmon,lit(js,'Std2Exponents')));E8=list(map(evalmon,lit(js,'Std8Exponents')));assert len(E2)==75 and len(E8)==375 and rank(E2)==75 and rank(E8)==375
Delta=sum(1<<(i*i)for i in range(1,25,2)if i*i<N);DE2=[mul(Delta,x)for x in E2]
assert len(E2)+2*len(E8)-rank([a|(b<<600)for a,b in zip(E2,DE2)]+E8+[x<<600 for x in E8])==17
G=combine(lit(js,'Gvec'),E2);A=combine(lit(js,'NumeratorVec'),E8);B=combine(lit(js,'DenominatorVec'),E8);assert G and A==G and B==mul(Delta,G)and B
print('NONIC71_Q_IDENTITIES_READONLY_PASS canonical_degree8 auxiliary_degree2 simultaneous_dimension17 q_precision600 Sturm480 quadrics276 diamond7_order5')
