#!/usr/bin/env python3
"""Independent canonical R2 multiplication replay for64 C87 degree9 cases."""
from pathlib import Path
import ast,re,json,gzip,hashlib,sys,time,random
from sage.all import GF,matrix
H=Path(__file__).resolve().parents[1]/'data'/'canonical_replay';C=H.parent;Q=C/'quotient33';t0=time.monotonic();REPLAY='--readonly'in sys.argv

def lit(s,n):return ast.literal_eval(re.search(re.escape(n)+r':=(.*?);',s,re.S)[1].replace('[*','[').replace('*]',']').replace('<','(').replace('>',')'))
def matbits(s,n):
 m=re.search(re.escape(n)+r':=Matrix\((?:GF\(2\)|Integers\(\)),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,n;r,c=map(int,m.group(1,2));a=ast.literal_eval(m[3]);assert len(a)==r*c;return[sum((int(a[i*c+j])&1)<<j for j in range(c))for i in range(r)]
def bits(a):return sum(int(v)<<i for i,v in enumerate(a))
def xorrows(v,rows):
 z=0
 while v:k=v&-v;z^=rows[k.bit_length()-1];v^=k
 return z
class Span:
 def __init__(self,rows=()):
  self.p={}
  for a in rows:self.add(a)
 def add(self,a):
  while a:
   k=(a&-a).bit_length()-1
   if k in self.p:a^=self.p[k]
   else:self.p[k]=a;return True
  return False
 def quotient(self,a):
  for k in sorted(self.p):
   if a>>k&1:a^=self.p[k]
  return a
 def __len__(self):return len(self.p)
def kernel(rows):
 p={};out=[]
 for i,a in enumerate(rows):
  c=1<<i
  while a:
   k=(a&-a).bit_length()-1
   if k in p:b,v=p[k];a^=b;c^=v
   else:p[k]=(a,c);break
  else:out.append(c)
 return out
class Field:
 def __init__(self,p):self.p=p;self.d=p.bit_length()-1
 def mul(self,a,b):
  z=0
  while b:
   if b&1:z^=a
   a<<=1;b>>=1
  while z.bit_length()>self.d:z^=self.p<<(z.bit_length()-1-self.d)
  return z
 def pow(self,a,n):
  z=1
  while n:
   if n&1:z=self.mul(z,a)
   a=self.mul(a,a);n>>=1
  return z
 def inv(self,a):assert a;return self.pow(a,(1<<self.d)-2)
 def smul(self,a,b,n=6):
  c=[0]*n
  for i,x in enumerate(a[:n]):
   for j,y in enumerate(b[:n-i]):
    if x and y:c[i+j]^=self.mul(x,y)
  return c
 def inverse_matrix(self,A):
  n=len(A);M=[row[:]+[int(i==j)for j in range(n)]for i,row in enumerate(A)]
  for k in range(n):
   r=next(i for i in range(k,n)if M[i][k]);M[k],M[r]=M[r],M[k];iv=self.inv(M[k][k]);M[k]=[self.mul(x,iv)for x in M[k]]
   for i in range(n):
    if i!=k and M[i][k]:a=M[i][k];M[i]=[x^self.mul(a,y)for x,y in zip(M[i],M[k])]
  assert[row[:n]for row in M]==[[int(i==j)for j in range(n)]for i in range(n)];return[row[n:]for row in M]
 def apply(self,A,v):
  out=[]
  for row in A:
   a=0
   for x,y in zip(row,v):a^=self.mul(x,y)
   out.append(a)
  return out
cs=(Q/'canonical_data.m').read_text();js=(Q/'jmap_q_data.m').read_text();sp=(C/'source_places_data.m').read_text();tp=(C/'target_places_data.m').read_text();q=matbits(js,'qcoeffs900');relations=matbits(cs,'QuadRelations');pairs=lit(cs,'pairs');assert len(q)==33 and len(relations)==465 and len(pairs)==561;assert len(Span(q))==33;assert [v&((1<<190)-1)for v in q]==matbits(cs,'qcoeffs');mask=(1<<900)-1

def qmul(a,b):
 if a.bit_count()<b.bit_count():a,b=b,a
 z=0
 while b:k=b&-b;z^=a<<(k.bit_length()-1);b^=k
 return z&mask
qq=[qmul(q[i-1],q[j-1])for i,j in pairs];assert len(Span(qq))==96
for rel in relations:assert xorrows(rel,qq)==0
# Select a basis of R2 and express all R2xR2 products in R4 by q expansions.
piv={};e2=[];rq=[]
for pair,v in zip(pairs,qq):
 w=v
 while w:
  k=w.bit_length()-1
  if k in piv:w^=piv[k]
  else:piv[k]=w;e2.append(pair);rq.append(v);break
assert len(e2)==96
piv={};r4=[];TT=[[0]*96 for _ in range(96)]
for i in range(96):
 for j in range(i,96):
  v=qmul(rq[i],rq[j]);w=v;c=0
  while w:
   k=w.bit_length()-1
   if k in piv:a,b=piv[k];w^=a;c^=b
   else:n=len(r4);piv[k]=(w,c^(1<<n));r4.append(v);c=1<<n;break
  TT[i][j]=TT[j][i]=c
assert len(r4)==224
for i in range(96):
 for j in range(96):assert xorrows(TT[i][j],r4)==qmul(rq[i],rq[j])
# Every retained finite point satisfies all465 canonical quadrics.
sk=lit(sp,'SourcePlaceKeys');tk=lit(tp,'TargetPlaceKeys');reps=lit(sp,'SourceRepresentativeIndices');assert len(sk)==len(tk)==4 and all(len(v)==4 for v in sk)and all(len(v)==16 for v in tk)and reps==[1]*4

def check_key(key):
 p,coords=key;F=Field(bits(p));x=[bits(a)for a in coords];assert len(x)==33 and next(a for a in x if a)==1;vals=[F.mul(x[i-1],x[j-1])for i,j in pairs]
 for rel in relations:assert xorrows(rel,vals)==0
 # Exact irreducibility follows from the finite-field Frobenius criterion.
 a=2
 for _ in range(F.d):a=F.mul(a,a)
 assert a==2
 for prime in [n for n in range(2,F.d+1)if F.d%n==0 and all(n%d for d in range(2,int(n**.5)+1))]:
  a=2
  for _ in range(F.d//prime):a=F.mul(a,a)
  b=a^2;c=F.p
  while b:
   aa=c
   while aa.bit_length()>=b.bit_length():aa^=b<<(aa.bit_length()-b.bit_length())
   c,b=b,aa
  assert c==1
 return F,x
for group in sk+tk:
 for key in group:check_key(key)
JETS=[];metadata=[]
for sid,group in enumerate(sk,1):
 F,x=check_key(group[0]);assert F.d==9;k=next(i for i,a in enumerate(x)if a);grads=[]
 for rel in relations:
  row=[0]*33
  while rel:
   b=rel&-rel;i,j=pairs[b.bit_length()-1];rel^=b
   if i!=j:row[i-1]^=x[j-1];row[j-1]^=x[i-1]
  grads.append(row)
 cols=[i for i in range(33)if i!=k];piv={};selected=[]
 for ix,g in enumerate(grads):
  row=[g[i]for i in cols]
  for j,p in sorted(piv.items()):
   if row[j]:a=row[j];row=[u^F.mul(a,v)for u,v in zip(row,p)]
  nz=next((j for j,a in enumerate(row)if a),None)
  if nz is not None:iv=F.inv(row[nz]);piv[nz]=[F.mul(a,iv)for a in row];selected.append(ix)
 assert len(piv)==31;free=next(i for i in range(32)if i not in piv);l=cols[free];unknown=[i for i in cols if i!=l];inv=F.inverse_matrix([[grads[ix][j]for j in unknown]for ix in selected]);jets=[[a]+[0]*5 for a in x]
 for n in range(1,6):
  jets[l][n]=int(n==1);vals=[]
  for i,j in pairs:
   v=0
   for a in range(n+1):v^=F.mul(jets[i-1][a],jets[j-1][n-a])
   vals.append(v)
  sol=F.apply(inv,[xorrows(relations[ix],vals)for ix in selected])
  for j,a in zip(unknown,sol):jets[j][n]=a
  vals=[]
  for i,j in pairs:
   v=0
   for a in range(n+1):v^=F.mul(jets[i-1][a],jets[j-1][n-a])
   vals.append(v)
  assert all(xorrows(rel,vals)==0 for rel in relations)
 assert any(a[1]for a in jets);JETS.append((F,jets));metadata.append({'source':sid,'field_polynomial':F.p,'normalized_coordinate':k+1,'uniformizer_coordinate':l+1,'independent_quadric_rows':[i+1 for i in selected],'Jacobian_rank':31,'precision':6,'jets':jets});print('C87_CANONICAL_JET_REPLAY source',sid,'rank31 precision6 all465quadrics=true',flush=True)
# Riemann--Roch injection, with B=Q54-6P9.
rng=random.Random(20261026);saved=json.loads(gzip.open(H/'matrices64.json.gz','rt').read())if REPLAY else None;MATRICES=[];RESULTS=[];UNRES=[];ix=0
for sid,(F,jets)in enumerate(JETS,1):
 bm=[sum(c<<(9*n)for n,c in enumerate(F.smul(jets[i-1],jets[j-1])))for i,j in e2]
 for qi,key in enumerate(tk[sid-1],1):
  FF,z=check_key(key);assert FF.d==54;bp=[FF.mul(z[i-1],z[j-1])for i,j in e2];kk=kernel(bp);assert len(kk)>=42;old=saved[ix]if REPLAY else None;fv=int(old['F'],16)if REPLAY else kk[0];assert fv and xorrows(fv,bp)==0
  W=Span([xorrows(fv,[row[j]for row in TT])for j in range(96)]);assert len(W)==96;free=[k for k in range(224)if k not in W.p];assert len(free)==128;rows=bm[:];hh=[];cand=[int(a,16)for a in old['H']]if REPLAY else[xorrows(rng.getrandbits(len(kk)),kk)for _ in range(16)]
  for h in cand:
   if not h or h==fv:continue
   assert xorrows(h,bp)==0;hh.append(hex(h))
   for i in range(96):w=W.quotient(xorrows(h,TT[i]));rows[i]|=sum(((w>>j)&1)<<k for k,j in enumerate(free))<<(54+128*(len(hh)-1))
   rank=len(Span(rows))
   if rank==96:break
  rec={'source':sid,'target':qi,'F':hex(fv),'H':hh,'columns':54+128*len(hh),'rank':rank}
  if rank!=96:assert not REPLAY;UNRES.append(rec)
  else:
   mm={k:rec[k]for k in['source','target','F','H','columns']}|{'rows':[hex(a)for a in rows]}
   if REPLAY:assert mm==old
   assert matrix(GF(2),[[a>>j&1 for j in range(rec['columns'])]for a in rows]).rank()==96;MATRICES.append(mm);RESULTS.append(rec)
  ix+=1
assert ix==64
inputs={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()for p in[Q/'canonical_data.m',Q/'jmap_q_data.m',C/'source_places_data.m',C/'target_places_data.m']};out={'status':'PARTIAL'if UNRES else'PASS','cases':len(RESULTS),'unresolved':UNRES,'genus':33,'R2_dimension':96,'R4_dimension':224,'source_formal_jet_precision':6,'source_Jacobian_rank':31,'quadrics':465,'products':9216,'q_precision':900,'Sturm_weight8':320,'inputs':inputs,'results':RESULTS}
if REPLAY:old_report=json.loads((H/'verified64.json').read_text());assert sorted(old_report.pop('inputs').values())==sorted(out['inputs'].values())and old_report=={k:v for k,v in out.items()if k!='inputs'};assert len(saved)==64
else:
 (H/'verified64.json').write_text(json.dumps(out,indent=2)+'\n');(H/'source_formal_jets.json').write_text(json.dumps(metadata,indent=2)+'\n')
 with gzip.open(H/'matrices64.json.gz','wt')as f:json.dump(MATRICES,f,separators=(',',':'))
print('C87_CANONICAL64_'+('READONLY_'if REPLAY else'')+'RESULT',len(RESULTS),'excluded',len(UNRES),'unresolved elapsed',time.monotonic()-t0,flush=True)
