#!/usr/bin/env python3
"""Exact canonical-ring products, certified by q-expansion beyond Sturm."""
from pathlib import Path
import ast,re,json,time,hashlib
H=Path(__file__).resolve().parents[1]/'data'/'quotient19';t0=time.monotonic();s=(H/'jmap_q_data.m').read_text();N=600;MASK=(1<<N)-1

def lit(n):return ast.literal_eval(re.search(re.escape(n)+r':=(.*?);',s,re.S)[1].replace('<','(').replace('>',')'))
def matrix(n):
 m=re.search(re.escape(n)+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m
 r,c=map(int,m.group(1,2));a=ast.literal_eval(m[3]);assert len(a)==r*c
 return [sum(a[i*c+j]<<j for j in range(c))for i in range(r)]
def mul(a,b):
 if a.bit_count()<b.bit_count():a,b=b,a
 z=0
 while b:
  k=b&-b;z^=a<<(k.bit_length()-1);b^=k
 return z&MASK

def insert(piv,x,comb=0):
 while x:
  k=(x&-x).bit_length()-1
  if k in piv:a,b=piv[k];x^=a;comb^=b
  else:piv[k]=(x,comb);return True
 return False

def pivots(rows):
 p={}
 for i,x in enumerate(rows):assert insert(p,x,1<<i)
 return p

def coordinates(x,p):
 z=0
 while x:
  k=(x&-x).bit_length()-1;assert k in p
  a,b=p[k];x^=a;z^=b
 return z


q=matrix('qcoeffs600');assert len(q)==19
powers=[[1]for _ in q]
for i,f in enumerate(q):
 for k in range(1,7):powers[i].append(mul(powers[i][-1],f))
def mon(e):
 x=1
 for i,k in enumerate(e):
  if k:x=mul(x,powers[i][k])
 return x
E3=lit('Std3Exponents');Q3=list(map(mon,E3));assert len(E3)==90;P3=pivots(Q3)
E6=[];Q6=[];p={};seen=set();done=False
for i,a in enumerate(E3):
 for j,b in enumerate(E3):
  e=tuple(x+y for x,y in zip(a,b))
  if e in seen:continue
  seen.add(e);w=mul(Q3[i],Q3[j])
  if insert(p,w):E6.append(list(e));Q6.append(w)
  if len(E6)==198:done=True;break
 if done:break
assert done;P6=pivots(Q6);TT=[]
for i,a in enumerate(Q3):
 TT.append([hex(coordinates(mul(a,b),P6))for b in Q3])
 if i%20==0:print('NONIC73_CANONICAL_TENSOR_PROGRESS',i,time.monotonic()-t0,flush=True)
stored=json.loads((H/'canonical_ring_tensor.json').read_text());assert stored['degree3_exponents']==E3 and stored['degree6_exponents']==E6 and stored['product_3_3_into_6']==TT
print('NONIC73_CANONICAL_RING_TENSOR_READONLY_PASS products8100 q_precision600 Sturm296')
