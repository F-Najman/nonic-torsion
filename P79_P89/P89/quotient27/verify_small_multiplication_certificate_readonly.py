#!/usr/bin/env python3
"""Exact canonical-ring Riemann--Roch sieve using a residual Artin algebra."""
from pathlib import Path
import ast,re,json,sys,time,hashlib
H=Path(__file__).resolve().parents[2]/'data/P89/quotient27';t0=time.monotonic();MATRICES=[]
def lit(s,n):return ast.literal_eval(re.search(re.escape(n)+r':=(.*?);',s,re.S)[1].replace('[*','[').replace('*]',']').replace('<','(').replace('>',')'))
def mat(s,n):
 m=re.search(re.escape(n)+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,n
 r,c=map(int,m.group(1,2));a=ast.literal_eval(m[3]);assert len(a)==r*c
 return [sum(a[i*c+j]<<j for j in range(c))for i in range(r)]
def bits(a):return sum(x<<i for i,x in enumerate(a))
def xor_rows(v,rows):
 z=0
 while v:
  b=v&-v;z^=rows[b.bit_length()-1];v^=b
 return z
class Span:
 def __init__(self,rows=()):
  self.p={}
  for x in rows:self.add(x)
 def reduce(self,x):
  while x:
   k=(x&-x).bit_length()-1
   if k not in self.p:break
   x^=self.p[k]
  return x
 def add(self,x):
  x=self.reduce(x)
  if not x:return False
  self.p[(x&-x).bit_length()-1]=x;return True
 # Fully reduce in all pivot columns, retaining nonpivot columns.
 def quotient(self,x):
  for k in sorted(self.p):
   if x>>k&1:x^=self.p[k]
  return x
 def __len__(self):return len(self.p)
def kernel(rows):
 p={};out=[]
 for i,x in enumerate(rows):
  c=1<<i
  while x:
   k=(x&-x).bit_length()-1
   if k in p:a,b=p[k];x^=a;c^=b
   else:p[k]=(x,c);break
  else:out.append(c)
 return out
class Field:
 def __init__(self,p):self.p=p;self.d=p.bit_length()-1
 def mul(self,a,b):
  z=0
  while b:
   if b&1:z^=a
   b>>=1;a<<=1
  while z.bit_length()>self.d:z^=self.p<<(z.bit_length()-1-self.d)
  return z
 def pow(self,a,n):
  z=1
  while n:
   if n&1:z=self.mul(z,a)
   a=self.mul(a,a);n>>=1
  return z
 def inv(self,a):assert a;return self.pow(a,(1<<self.d)-2)
 def smul(self,a,b,n):
  c=[0]*n
  for i,x in enumerate(a):
   for j,y in enumerate(b[:n-i]):
    if x and y:c[i+j]^=self.mul(x,y)
  return c
 def sinv(self,a):
  n=len(a);iv=self.inv(a[0]);b=[iv]
  for k in range(1,n):
   x=0
   for i in range(1,k+1):x^=self.mul(a[i],b[k-i])
   b.append(self.mul(iv,x))
  return b
 def normalize(self,z):
  k=next(i for i,x in enumerate(z)if x[0]);iv=self.sinv(z[k]);return [self.smul(x,iv,len(iv))for x in z]
 def rank(self,rows):
  p={}
  for row in rows:
   row=list(row)
   for k,a in sorted(p.items()):
    if row[k]:
     c=row[k];row=[x^self.mul(c,y)for x,y in zip(row,a)]
   nz=next((i for i,x in enumerate(row)if x),None)
   if nz is not None:
    iv=self.inv(row[nz]);p[nz]=[self.mul(x,iv)for x in row]
  return len(p)
T=json.loads((H/'canonical_ring_tensor.json').read_text());T32=[[int(v,16)for v in row]for row in T['product_3_2_into_5']];T35=[[int(v,16)for v in row]for row in T['product_3_5_into_8']];E3=T['degree3_exponents'];E5=T['degree5_exponents'];assert T['dimensions']==[78,130,234,390]
ss=(H/'source_places_data.m').read_text();src_keys=lit(ss,'SourcePlaceKeys');perm=lit(ss,'SourceSigmaPermutation');orbs=lit(ss,'SourceSigmaOrbits');ji=lit(ss,'SourceJIndices')
cs=(H/'canonical_data.m').read_text();relations=mat(cs,'QuadRelations');pairs=lit(cs,'pairs');D=mat((H/'jmap_q_data.m').read_text(),'Diamond3Matrix')
jetrecords=lit((H/'canonical_formal_source_jets_v2_data.m').read_text(),'CanonicalSourceJets');jets={};fields={}
for rep,pol,coords in jetrecords:
 F=Field(bits(pol));z=[[bits(a)for a in c]for c in coords];assert F.d==9 and len(z)==27 and all(len(a)==5 for a in z)
 assert [a[0]for a in z]==[bits(a)for a in src_keys[rep-1][1]];assert any(a[1]for a in z)
 products=[F.smul(z[i-1],z[j-1],5)for i,j in pairs]
 for rel in relations:
  out=[0]*5
  while rel:
   b=rel&-rel;a=products[b.bit_length()-1];out=[x^y for x,y in zip(out,a)];rel^=b
  assert not any(out)
 grads=[]
 for rel in relations:
  row=[0]*27
  while rel:
   b=rel&-rel;i,j=pairs[b.bit_length()-1];rel^=b
   if i!=j:row[i-1]^=z[j-1][0];row[j-1]^=z[i-1][0]
  grads.append(row)
 assert F.rank(grads)==25
 k=rep;original=z
 for _ in range(4):
  assert k not in jets;jets[k]=z;fields[k]=F
  assert [a[0]for a in z]==[bits(a)for a in src_keys[k-1][1]]
  sz=[]
  for row in D:
   zz=[0]*5
   while row:
    b=row&-row;zz=[x^y for x,y in zip(zz,z[b.bit_length()-1])];row^=b
   sz.append(zz)
  z=F.normalize(sz);k=perm[k-1]
 assert k==rep and z==original
assert len(jets)==24
print('NONIC89_CANONICAL_JETS_INDEPENDENT_PASS',len(jetrecords),'smooth_points_and_transported_jets',len(jets),flush=True)

def eval_basis(z,F,exps,n):
 pp=[]
 for a in z:
  powers=[[1]+[0]*(n-1)]
  for k in range(1,6):powers.append(F.smul(powers[-1],a[:n],n))
  pp.append(powers)
 rows=[]
 for es in exps:
  w=[1]+[0]*(n-1)
  for i,k in enumerate(es):
   if k:w=F.smul(w,pp[i][k],n)
  rows.append(sum(a<<(F.d*k)for k,a in enumerate(w)))
 return rows
sev={}
for k,z in jets.items():sev[k]=[eval_basis(z,fields[k],E3,5),eval_basis(z,fields[k],E5,5)]
tev={};tperm={}
for f in(1,2,3):
 s=(H/f'target_places{f}_jv_data.m').read_text();ks=lit(s,'TargetPlaceKeys');tperm[f]=lit(s,'TargetSigmaPermutation')
 for k,key in enumerate(ks,1):
  F=Field(bits(key[0]));z=[[bits(c)]for c in key[1]];assert F.d==54
  tev[f,k]=[eval_basis(z,F,E3,1),eval_basis(z,F,E5,1)]
print('NONIC89_CANONICAL_LOCAL_EVALUATIONS_READY',time.monotonic()-t0,flush=True)

def block(f,target,srcsingle,srcfive,which):
 a=tev[f,target][which];b=sev[srcsingle][which];c=sev[srcfive][which]
 return [x|((y&511)<<54)|(z<<63)for x,y,z in zip(a,b,c)]

import gzip
with gzip.open(H/'canonical_ring_small_matrices.json.gz','rt')as fp:cases=json.load(fp)
expected=json.loads((H/'canonical_ring_small_sieve_all_verified.json').read_text());assert expected['cases']==96 and len(cases)==96
coverage={(ji[o[0]-1],o[0],q,c)for o in orbs for q in range(1,9)for c in(1,2)};seen=set()
for data,rec in zip(cases,expected['results']):
 f,p,q,convention=[data[k]for k in ('fiber','source','target','convention')];case=(f,p,q,convention);assert case in coverage and case not in seen;seen.add(case)
 p1=perm[p-1];p2=perm[p1-1];p3=perm[p2-1];q1=tperm[f][q-1]
 if convention==1:bp=[block(f,q1,p2,p,w)for w in(0,1)];bm=block(f,q,p3,p1,0)
 else:bp=[block(f,q1,p,p2,w)for w in(0,1)];bm=block(f,q,p1,p3,0)
 F,K=int(data['F'],16),int(data['K'],16);assert F and F.bit_length()<=130 and K.bit_length()<=234
 assert xor_rows(F,bp[0])==0 and xor_rows(K,bp[1])==0
 W=Span([xor_rows(F,[row[j]for row in T35])for j in range(234)]);assert len(W)==234
 free=[k for k in range(390)if k not in W.p];assert len(free)==156
 rows=[]
 for i in range(130):
  z=W.quotient(xor_rows(K,T35[i]));packed=sum(((z>>p)&1)<<h for h,p in enumerate(free));rows.append(bm[i]|(packed<<108))
 assert data['columns']==264 and data['rows']==[hex(v)for v in rows]
 assert hashlib.sha256('\n'.join(data['rows']).encode()).hexdigest()==rec['matrix_sha256']
 piv={}
 for x in rows:
  while x:
   k=x.bit_length()-1
   if k in piv:x^=piv[k]
   else:piv[k]=x;break
 assert len(piv)==130
assert seen==coverage
print('NONIC89_SMALL_MULTIPLICATION_CERTIFICATE_READONLY_PASS cases=96 rank=130 matrix=130x264 explicit_sections_vanish_positive=true time=',time.monotonic()-t0)
