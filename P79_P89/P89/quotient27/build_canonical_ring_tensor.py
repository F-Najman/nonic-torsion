#!/usr/bin/env python3
"""Exact canonical-ring products, certified by q-expansion beyond Sturm."""
from pathlib import Path
import ast,re,json,time,hashlib
H=Path(__file__).resolve().parents[2]/'data/P89/quotient27';t0=time.monotonic();s=(H/'jmap_q_data.m').read_text();N=600;MASK=(1<<N)-1

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
q=matrix('qcoeffs600');assert len(q)==27
powers=[[1]for _ in q]
for i,f in enumerate(q):
 for k in range(1,9):powers[i].append(mul(powers[i][-1],f))
def mon(e):
 x=1
 for i,k in enumerate(e):
  if k:x=mul(x,powers[i][k])
 return x
E2=lit('Std2Exponents');E8=lit('Std8Exponents');Q2=list(map(mon,E2));Q8=list(map(mon,E8));assert len(E2)==78 and len(E8)==390
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
E1=[[int(j==i)for j in range(27)]for i in range(27)]
E3,Q3,P3=basis_products(E1,q,E2,Q2,130);E5,Q5,P5=basis_products(E3,Q3,E2,Q2,234)
print('NONIC89_CANONICAL_RING_BASES_READY',len(E3),len(E5),flush=True)
T32=[];T35=[]
for i,a in enumerate(Q3):
 T32.append([hex(coordinates(mul(a,b),P5))for b in Q2])
 T35.append([hex(coordinates(mul(a,b),P8))for b in Q5])
 if i%20==0:print('NONIC89_CANONICAL_RING_TENSOR_PROGRESS',i,'time',time.monotonic()-t0,flush=True)
record={'genus':27,'q_precision':600,'weight16_Sturm_bound':480,'degrees':[2,3,5,8],'dimensions':[78,130,234,390],'degree2_exponents':E2,'degree3_exponents':E3,'degree5_exponents':E5,'degree8_exponents':E8,'product_3_2_into_5':T32,'product_3_5_into_8':T35,'q_input_sha256':hashlib.sha256(s.encode()).hexdigest()}
(H/'canonical_ring_tensor.json').write_text(json.dumps(record,separators=(',',':'))+'\n')
summary={'status':'PASS','all_q_product_identities_verified':True,'product_3_2_count':130*78,'product_3_5_count':130*234,'q_precision':600,'Sturm_bound':480,'dimensions':[78,130,234,390],'elapsed_seconds':time.monotonic()-t0,'output_sha256':hashlib.sha256((H/'canonical_ring_tensor.json').read_bytes()).hexdigest()}
(H/'canonical_ring_tensor_verified.json').write_text(json.dumps(summary,indent=2)+'\n');print('NONIC89_CANONICAL_RING_TENSOR_COMPLETE',json.dumps(summary,sort_keys=True))
