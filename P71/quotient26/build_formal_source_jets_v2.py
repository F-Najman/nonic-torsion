#!/usr/bin/env python3
"""Exact formal parameters from smooth canonical quadrics overF512."""
from pathlib import Path
import re,ast,json,hashlib,time
H=Path(__file__).resolve().parents[1]/'data'/'quotient26';t0=time.monotonic()
def lit(s,n):return ast.literal_eval(re.search(re.escape(n)+r':=(.*?);',s,re.S)[1].replace('[*','[').replace('*]',']').replace('<','(').replace('>',')'))
def bits(a):return sum(x<<i for i,x in enumerate(a))
def mat(s,n):
 m=re.search(re.escape(n)+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m
 r,c=map(int,m.group(1,2));a=ast.literal_eval(m[3]);return [sum(a[i*c+j]<<j for j in range(c))for i in range(r)]
class Field:
 def __init__(self,p):self.p=p;self.d=p.bit_length()-1
 def mul(self,a,b):
  z=0
  while b:
   if b&1:z^=a
   a<<=1;b>>=1
  while z.bit_length()>self.d:z^=self.p<<(z.bit_length()-1-self.d)
  return z
 def inv(self,a):
  assert a;n=(1<<self.d)-2;b=1
  while n:
   if n&1:b=self.mul(b,a)
   a=self.mul(a,a);n>>=1
  return b
 def inverse_matrix(self,A):
  n=len(A);M=[row[:]+[int(i==j)for j in range(n)]for i,row in enumerate(A)]
  for k in range(n):
   r=next(i for i in range(k,n)if M[i][k]);M[k],M[r]=M[r],M[k];iv=self.inv(M[k][k]);M[k]=[self.mul(x,iv)for x in M[k]]
   for i in range(n):
    if i!=k and M[i][k]:
     a=M[i][k];M[i]=[x^self.mul(a,y)for x,y in zip(M[i],M[k])]
  assert [row[:n]for row in M]==[[int(i==j)for j in range(n)]for i in range(n)]
  return [row[n:]for row in M]
 def apply(self,A,v):
  out=[]
  for row in A:
   z=0
   for a,b in zip(row,v):z^=self.mul(a,b)
   out.append(z)
  return out
s=(H/'source_places_data.m').read_text();keys=lit(s,'SourcePlaceKeys');orbs=lit(s,'SourceSigmaOrbits');c=(H/'canonical_data.m').read_text();rels=mat(c,'QuadRelations');pairs=lit(c,'pairs');records=[];metadata=[]
def xor(v,values):
 z=0
 while v:
  k=v&-v;z^=values[k.bit_length()-1];v^=k
 return z
for orb in orbs:
 rep=orb[0];key=keys[rep-1];F=Field(bits(key[0]));assert F.d in[7,8];x=[bits(a)for a in key[1]];k=next(i for i,a in enumerate(x)if a);assert x[k]==1
 grads=[]
 for rel in rels:
  row=[0]*26
  while rel:
   b=rel&-rel;i,j=pairs[b.bit_length()-1];rel^=b
   if i!=j:row[i-1]^=x[j-1];row[j-1]^=x[i-1]
  grads.append(row)
 cols=[i for i in range(26)if i!=k];piv={};selected=[]
 for ix,g in enumerate(grads):
  row=[g[i]for i in cols]
  for j,p in sorted(piv.items()):
   if row[j]:a=row[j];row=[u^F.mul(a,v)for u,v in zip(row,p)]
  nz=next((j for j,a in enumerate(row)if a),None)
  if nz is not None:
   iv=F.inv(row[nz]);piv[nz]=[F.mul(a,iv)for a in row];selected.append(ix)
 assert len(piv)==24 and len(selected)==24
 free=next(i for i in range(25)if i not in piv);l=cols[free];unknown=[i for i in cols if i!=l];inv=F.inverse_matrix([[grads[ix][j]for j in unknown]for ix in selected]);jets=[[a,0,0,0,0]for a in x]
 for n in range(1,5):
  jets[l][n]=int(n==1)
  products=[]
  for i,j in pairs:
   z=0
   for a in range(n+1):z^=F.mul(jets[i-1][a],jets[j-1][n-a])
   products.append(z)
  b=[xor(rels[ix],products)for ix in selected];solution=F.apply(inv,b)
  for j,a in zip(unknown,solution):jets[j][n]=a
  products=[]
  for i,j in pairs:
   z=0
   for a in range(n+1):z^=F.mul(jets[i-1][a],jets[j-1][n-a])
   products.append(z)
  assert all(xor(rel,products)==0 for rel in rels)
 assert jets[k]==[1,0,0,0,0]and jets[l][1:]==[1,0,0,0]
 encoded=[[[a>>r&1 for r in range(F.d)]for a in row]for row in jets];records.append((rep,key[0],encoded));metadata.append({'source_representative':rep,'normalized_coordinate':k+1,'uniformizer_coordinate':l+1,'independent_equations':[i+1 for i in selected],'Jacobian_rank':24,'jet_precision':5})
 print('NONIC71_FORMAL_CANONICAL_JETS_POINT',rep,'smooth_rank24','uniformizer_coordinate',l+1,flush=True)
body='CanonicalSourceJets:=[* '+',\n'.join('<'+str(r)+','+str(p)+','+str(a)+'>'for r,p,a in records)+' *];CanonicalJetPrecision:=5;\nprint "NONIC71_CANONICAL_FORMAL_SOURCE_JETS_READY";\n'
out=H/'canonical_formal_source_jets_v2_data.m';out.write_text(body)
record={'status':'PASS','representatives':2,'all_quadric_jet_equations_verified':276,'precision':5,'records':metadata,'canonical_data_sha256':hashlib.sha256(c.encode()).hexdigest(),'source_places_sha256':hashlib.sha256(s.encode()).hexdigest(),'jet_data_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t0}
(H/'canonical_formal_source_jets_v2_verified.json').write_text(json.dumps(record,indent=2)+'\n');print('NONIC71_FORMAL_CANONICAL_JETS_COMPLETE',time.monotonic()-t0)
