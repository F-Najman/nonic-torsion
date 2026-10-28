#!/usr/bin/env python3
"""Exact canonical-ring Riemann--Roch sieve using a residual Artin algebra."""
from pathlib import Path
import ast,re,json,sys,time,hashlib
H=Path(__file__).resolve().parents[1]/'data'/'quotient45';t0=time.monotonic();MATRICES=[]
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


import random,gzip,itertools
REPLAY='--readonly'in sys.argv;rng=random.Random(20261026)
T=json.loads((H/'canonical_ring_tensor.json').read_text());TT=[[int(v,16)for v in row]for row in T['product_2_2_into_4']];E4=T['degree2_exponents'];assert T['dimensions']==[132,308]
assert lit((H/'finite_image_winding_data.m').read_text(),'AlphaExponent')==1
ss=(H/'source_places_data.m').read_text();src_keys=lit(ss,'SourcePlaceKeys');perm=lit(ss,'SourceSigmaPermutation');orbs=lit(ss,'SourceSigmaOrbits');ji=lit(ss,'SourceJIndices');ds=lit(ss,'SourcePlaceDegrees')
cs=(H/'canonical_data.m').read_text();relations=mat(cs,'QuadRelations');pairs=lit(cs,'pairs');D=mat((H/'jmap_q_data.m').read_text(),'Diamond2Matrix');jetrecords=lit((H/'canonical_formal_source_jets_v2_data.m').read_text(),'CanonicalSourceJets');jets={};fields={}
for rep,pol,coords in jetrecords:
 F=Field(bits(pol));z=[[bits(a)for a in c]for c in coords];assert F.d==ds[rep-1]and len(z)==45 and all(len(a)==5 for a in z);assert[a[0]for a in z]==[bits(a)for a in src_keys[rep-1][1]]and any(a[1]for a in z)
 products=[F.smul(z[i-1],z[j-1],5)for i,j in pairs]
 for rel in relations:
  out=[0]*5
  while rel:
   b=rel&-rel;a=products[b.bit_length()-1];out=[x^y for x,y in zip(out,a)];rel^=b
  assert not any(out)
 grads=[]
 for rel in relations:
  row=[0]*45
  while rel:
   b=rel&-rel;i,j=pairs[b.bit_length()-1];rel^=b
   if i!=j:row[i-1]^=z[j-1][0];row[j-1]^=z[i-1][0]
  grads.append(row)
 assert F.rank(grads)==43
 k=rep;original=z
 for _ in range(11):
  assert k not in jets;jets[k]=z;fields[k]=F;assert[a[0]for a in z]==[bits(a)for a in src_keys[k-1][1]]
  sz=[]
  for row in D:
   zz=[0]*5
   while row:
    b=row&-row;zz=[x^y for x,y in zip(zz,z[b.bit_length()-1])];row^=b
   sz.append(zz)
  z=F.normalize(sz);k=perm[k-1]
 assert k==rep and z==original
assert len(jets)==55

def eval_basis(z,F,exps,n):
 pp=[]
 for a in z:
  powers=[[1]+[0]*(n-1)]
  for k in range(1,9):powers.append(F.smul(powers[-1],a[:n],n))
  pp.append(powers)
 rows=[]
 for es in exps:
  w=[1]+[0]*(n-1)
  for i,k in enumerate(es):
   if k:w=F.smul(w,pp[i][k],n)
  rows.append(sum(a<<(F.d*k)for k,a in enumerate(w)))
 return rows
sev={k:eval_basis(z,fields[k],E4,3)for k,z in jets.items()};td={};te={};tp={};tk={}
for fi in range(1,8):
 s=(H/f'target_places{fi}_jbase_data.m').read_text();ks=lit(s,'TargetPlaceKeys');tp[fi]=lit(s,'TargetSigmaPermutation');td[fi]=lit(s,'TargetDegree');te[fi]={};tk[fi]=ks
 for k,key in enumerate(ks,1):
  F=Field(bits(key[0]));assert F.d==td[fi];te[fi][k]=eval_basis([[bits(c)]for c in key[1]],F,E4,1)
source_targets=[ast.literal_eval(a)for a in re.findall(r'Append\(~SourceTargetIndices,(\[.*?\])\);',(H.parent/'quotient15'/'isogeny3_target_data.m').read_text())];assert len(source_targets)==5
def combine(parts):
 out=[0]*132;offset=0
 for rows,width in parts:
  assert len(rows)==132 and all(a.bit_length()<=width for a in rows)
  out=[a|(b<<offset)for a,b in zip(out,rows)];offset+=width
 return out,offset
saved=json.loads(gzip.open(H/'canonical_ring_matrices.json.gz','rt').read())if REPLAY else None
RESULTS=[];MATRICES=[];UNRESOLVED=[];caseindex=0
for orb in orbs:
 p=orb[0];pp=[p]
 for _ in range(10):pp.append(perm[pp[-1]-1])
 assert len(set(pp))==11 and perm[pp[-1]-1]==p
 p1,p6,p7=pp[1],pp[6],pp[7];d=ds[p-1];fis=source_targets[ji[p-1]-1];mask=(1<<d)-1
 for choices in itertools.product(range(1,12),repeat=len(fis)):
  targets=list(zip(fis,choices))
  # In the degree8 branch the degree8 target lies in the OTHER j-orbit.
  # Verify that no repeated support is hidden in concatenated jet conditions.
  for fi,q in targets:
   assert tk[fi][q-1]not in[sk for i,sk in enumerate(src_keys,1)if i in pp]
   assert tk[fi][tp[fi][q-1]-1]not in[sk for i,sk in enumerate(src_keys,1)if i in pp]
  targetpos=[(te[fi][tp[fi][q-1]],td[fi])for fi,q in targets];targetneg=[(te[fi][q],td[fi])for fi,q in targets]
  for conv in(1,2):
   if conv==1:
    bp,deg=combine(targetpos+[([a&mask for a in sev[p6]],d),(sev[p],3*d)])
    bm,dd=combine(targetneg+[([a&mask for a in sev[p7]],d),(sev[p1],3*d)])
   else:
    bp,deg=combine(targetpos+[([a&mask for a in sev[p]],d),(sev[p6],3*d)])
    bm,dd=combine(targetneg+[([a&mask for a in sev[p1]],d),(sev[p7],3*d)])
   assert deg==dd==8*d;kk=kernel(bp);assert len(kk)>=132-deg;old=saved[caseindex]if REPLAY else None
   if REPLAY:assert old['source']==p and old['targets']==[list(a)for a in targets]and old['convention']==conv
   fv=int(old['F'],16)if REPLAY else kk[0];assert fv and xor_rows(fv,bp)==0
   W=Span([xor_rows(fv,[row[j]for row in TT])for j in range(132)]);assert len(W)==132;free=[k for k in range(308)if k not in W.p];assert len(free)==176
   rows=bm[:];hh=[];candidates=[int(a,16)for a in old['H']]if REPLAY else[xor_rows(rng.getrandbits(len(kk)),kk)for _ in range(16)]
   for h in candidates:
    if not h or h==fv:continue
    assert xor_rows(h,bp)==0;hh.append(hex(h))
    for i in range(132):
     w=W.quotient(xor_rows(h,TT[i]));rows[i]|=sum(((w>>j)&1)<<k for k,j in enumerate(free))<<(deg+176*(len(hh)-1))
    rank=len(Span(rows))
    if rank==132:break
   rec={'source':p,'source_j_index':ji[p-1],'source_degree':d,'targets':[list(a)for a in targets],'convention':conv,'rank':rank,'F':hex(fv),'H':hh,'columns':deg+176*len(hh),'positive_degree':deg}
   if rank!=132:
    assert not REPLAY;UNRESOLVED.append(rec);caseindex+=1;continue
   mm={k:rec[k]for k in('source','targets','convention','F','H','columns')}|{'rows':[hex(a)for a in rows]}
   if REPLAY:assert mm==old
   RESULTS.append(rec);MATRICES.append(mm);caseindex+=1
assert caseindex==550
if REPLAY:assert len(saved)==550;print('NONIC67H3_MULTIPLICATION_CERTIFICATE_READONLY_PASS matrices550 rank132 all_reconstructed=true elapsed',time.monotonic()-t0)
else:
 out={'status':'PARTIAL'if UNRESOLVED else'PASS','cases':len(RESULTS),'unresolved':UNRESOLVED,'genus':45,'AF_degree':42,'finite_image_operator':'diamond2-1','canonical_source_jets':55,'smooth_source_rank':43,'q_precision':600,'product_Sturm_bound':498,'results':RESULTS,'elapsed_seconds':time.monotonic()-t0};(H/'canonical_ring_sieve_all_verified.json').write_text(json.dumps(out,indent=2)+'\n')
 with gzip.open(H/'canonical_ring_matrices.json.gz','wt')as fp:json.dump(MATRICES,fp,separators=(',',':'))
 print('NONIC67H3_CANONICAL_RING_SIEVE',len(RESULTS),'excluded',len(UNRESOLVED),'unresolved',time.monotonic()-t0,flush=True)
