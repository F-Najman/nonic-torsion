#!/usr/bin/env python3
"""Exact canonical-ring products, certified by q-expansion beyond Sturm."""
from pathlib import Path
import ast,re,json,time,hashlib
H=Path(__file__).resolve().parents[2]/'data/P79/quotient18';t0=time.monotonic();s=(H/'jmap_q_data.m').read_text();N=600;MASK=(1<<N)-1

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

q=matrix('qcoeffs600');assert len(q)==18
powers=[[1]for _ in q]
for i,f in enumerate(q):
 for k in range(1,9):powers[i].append(mul(powers[i][-1],f))
def mon(e):
 x=1
 for i,k in enumerate(e):
  if k:x=mul(x,powers[i][k])
 return x
E2=lit('Std2Exponents');E8=lit('Std8Exponents');Q2=list(map(mon,E2));Q8=list(map(mon,E8));assert len(E2)==51 and len(E8)==255
P2=pivots(Q2);P8=pivots(Q8)
def basis_products(A,QA,B,QB,target):
 ee=[];qq=[];p={};seen=set()
 for i,a in enumerate(A):
  for j,b in enumerate(B):
   e=tuple(x+y for x,y in zip(a,b))
   if e in seen:continue
   seen.add(e);w=mul(QA[i],QB[j])
   if insert(p,w):ee.append(list(e));qq.append(w)
   if len(ee)==target:return ee,qq,pivots(qq)
 raise AssertionError((len(ee),target))
E4,Q4,P4=basis_products(E2,Q2,E2,Q2,119);T44=[]
for i,a in enumerate(Q4):
 T44.append([hex(coordinates(mul(a,b),P8))for b in Q4])
 if i%20==0:print('NONIC79_CANONICAL_TENSOR_PROGRESS',i,time.monotonic()-t0,flush=True)
record={'genus':18,'q_precision':600,'weight16_Sturm_bound':320,'degrees':[4,8],'dimensions':[119,255],'degree4_exponents':E4,'degree8_exponents':E8,'product_4_4_into_8':T44,'q_input_sha256':hashlib.sha256(s.encode()).hexdigest()}
(H/'canonical_ring_tensor.json').write_text(json.dumps(record,separators=(',',':'))+'\n')
summary={'status':'PASS','all_q_product_identities_verified':True,'product_count':119*119,'q_precision':600,'Sturm_bound':320,'dimensions':[119,255],'elapsed_seconds':time.monotonic()-t0,'output_sha256':hashlib.sha256((H/'canonical_ring_tensor.json').read_bytes()).hexdigest()}
(H/'canonical_ring_tensor_verified.json').write_text(json.dumps(summary,indent=2)+'\n');print('NONIC79_CANONICAL_RING_TENSOR_COMPLETE',json.dumps(summary,sort_keys=True))
