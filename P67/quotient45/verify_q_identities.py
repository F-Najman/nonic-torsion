#!/usr/bin/env python3
"""Read-only exact q certificate for the degree9 canonical j-map."""
from pathlib import Path
import ast,re,json,hashlib
H=Path(__file__).resolve().parents[1]/'data'/'quotient45';cs=(H/'canonical_data.m').read_text();js=(H/'jmap_q_data.m').read_text();N=1400;mask=(1<<N)-1

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
old=matrix(cs,'qcoeffs');q=matrix(js,'qcoeffs1400');md=matrix(js,'Diamond2Matrix');rel=matrix(cs,'QuadRelations');assert old[:2]==(45,320)and q[:2]==(45,1400)and md[:2]==(45,45)and rel[:2]==(903,1035)
assert [x&((1<<320)-1)for x in q[2]]==old[2]and rank([x&((1<<125)-1)for x in old[2]])==45
I=[1<<i for i in range(45)];p=I
for n in range(1,12):
 p=[combine([(x>>j)&1 for j in range(45)],md[2])for x in p]
 if n<11:assert p!=I
assert p==I
pairs=lit(cs,'pairs');assert len(pairs)==1035 and set(pairs)=={(i,j)for i in range(1,46)for j in range(i,46)}
products=[mul(q[2][i-1],q[2][j-1])for i,j in pairs];assert rank(products)==132
for row in rel[2]:assert combine([(row>>j)&1 for j in range(1035)],products)==0
powers=[[1]for _ in range(45)]
for i in range(45):
 for k in range(1,11):powers[i].append(mul(powers[i][-1],q[2][i]))
def evalmon(e):
 assert len(e)==45;x=1
 for i,k in enumerate(e):
  if k:x=mul(x,powers[i][k])
 return x
E2=list(map(evalmon,lit(js,'Std4Exponents')));E8=list(map(evalmon,lit(js,'Std10Exponents')));assert len(E2)==308 and len(E8)==836 and rank(E2)==308 and rank(E8)==836
Delta=sum(1<<(i*i)for i in range(1,39,2)if i*i<N);DE2=[mul(Delta,x)for x in E2]
assert len(E2)+2*len(E8)-rank([a|(b<<1400)for a,b in zip(E2,DE2)]+E8+[x<<1400 for x in E8])==88
G=combine(lit(js,'Gvec'),E2);A=combine(lit(js,'NumeratorVec'),E8);B=combine(lit(js,'DenominatorVec'),E8);assert G and A==G and B==mul(Delta,G)and B
print('NONIC67H3_Q_IDENTITIES_READONLY_PASS canonical_degree10 auxiliary_degree4 simultaneous_dimension88 q_precision1400 Sturm1246 quadrics903 diamond2_order11')
