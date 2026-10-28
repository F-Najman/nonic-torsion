#!/usr/bin/env python3
"""Exact canonical-ring products, certified by q-expansion beyond Sturm."""
from pathlib import Path
import ast,re,json,time,hashlib
H=Path(__file__).resolve().parents[1]/'data'/'quotient26';t0=time.monotonic();s=(H/'jmap_q_data.m').read_text();N=600;MASK=(1<<N)-1

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


q=matrix('qcoeffs600');assert len(q)==26
powers=[[1]for _ in q]
for i,f in enumerate(q):
 for k in range(1,5):powers[i].append(mul(powers[i][-1],f))
def mon(e):
 x=1
 for i,k in enumerate(e):
  if k:x=mul(x,powers[i][k])
 return x
E2=lit('Std2Exponents');Q2=list(map(mon,E2));assert len(E2)==75;P2=pivots(Q2)
E4=[];Q4=[];p={};seen=set();done=False
for i,a in enumerate(E2):
 for j,b in enumerate(E2):
  e=tuple(x+y for x,y in zip(a,b))
  if e in seen:continue
  seen.add(e);w=mul(Q2[i],Q2[j])
  if insert(p,w):E4.append(list(e));Q4.append(w)
  if len(E4)==175:done=True;break
 if done:break
assert done;P4=pivots(Q4);TT=[]
for i,a in enumerate(Q2):
 TT.append([hex(coordinates(mul(a,b),P4))for b in Q2])
 if i%20==0:print('NONIC71_CANONICAL_TENSOR_PROGRESS',i,time.monotonic()-t0,flush=True)
stored=json.loads((H/'canonical_ring_tensor.json').read_text());assert stored['degree2_exponents']==E2 and stored['degree4_exponents']==E4 and stored['product_2_2_into_4']==TT
print('NONIC71_CANONICAL_RING_TENSOR_READONLY_PASS products5625 precision600 Sturm240')
