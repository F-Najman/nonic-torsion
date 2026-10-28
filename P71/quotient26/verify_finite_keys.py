#!/usr/bin/env python3
"""Read-only exact binary replay of canonical finite-place keys and diamonds."""
from pathlib import Path
import ast,re,json,sys
H=Path(__file__).resolve().parents[1]/'data'/'quotient26'

def lit(s,n):
 m=re.search(re.escape(n)+r':=(.*?);',s,re.S);assert m,n
 return ast.literal_eval(m[1].replace('[*','[').replace('*]',']').replace('<','(').replace('>',')').replace('true','True').replace('false','False'))
def mat(s,n):
 m=re.search(re.escape(n)+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,n
 r,c=map(int,m.group(1,2));a=ast.literal_eval(m[3]);assert len(a)==r*c
 return [a[i*c:(i+1)*c]for i in range(r)]
def bits(a):assert all(x in(0,1)for x in a);return sum(x<<i for i,x in enumerate(a))
def rem(a,p):
 d=p.bit_length()-1
 while a.bit_length()>d:a^=p<<(a.bit_length()-1-d)
 return a
def mul(a,b,p):
 z=0
 while b:
  if b&1:z^=a
  a<<=1;b>>=1
 return rem(z,p)
def powmod(a,n,p):
 z=1
 while n:
  if n&1:z=mul(z,a,p)
  a=mul(a,a,p);n>>=1
 return z
def gcd(a,b):
 while b:a,b=b,rem(a,b)
 return a
def irreducible(p):
 d=p.bit_length()-1;v=2
 for _ in range(d):v=mul(v,v,p)
 assert v==2
 for r in [q for q in range(2,d+1)if d%q==0 and all(q%a for a in range(2,int(q**.5)+1))]:
  v=2
  for _ in range(d//r):v=mul(v,v,p)
  assert gcd(v^2,p)==1

def norm(z,p):
 a=next(x for x in z if x);b=powmod(a,(1<<(p.bit_length()-1))-2,p)
 return [mul(x,b,p)for x in z]
def linear(row,z):
 x=0
 for a,b in zip(row,z):
  if a:x^=b
 return x
canon=(H/'canonical_data.m').read_text();qm=(H/'jmap_q_data.m').read_text();model=(H/'plane_model_certified_data.m').read_text()
D=mat(qm,'Diamond7Matrix');rels=mat(canon,'QuadRelations');pairs=lit(canon,'pairs');ee=lit(qm,'Std8Exponents');nv=lit(qm,'NumeratorVec');dv=lit(qm,'DenominatorVec');proj=mat(model,'ProjectionRows')
seen_polys=set()
def check(keys,perm=None,uv=None):
 zz=[];polys=[]
 for k,key in enumerate(keys):
  p=bits(key[0]);d=p.bit_length()-1
  if p not in seen_polys:irreducible(p);seen_polys.add(p)
  z=norm([bits(a)for a in key[1]],p);assert len(z)==26 and all(len(a)==d for a in key[1]);cz=[mul(z[i-1],z[j-1],p)for i,j in pairs]
  assert all(linear(row,cz)==0 for row in rels)
  powers=[[1]for _ in z]
  for i in range(26):
   for _ in range(8):powers[i].append(mul(powers[i][-1],z[i],p))
  ms=[]
  for es in ee:
   w=1
   for i,e in enumerate(es):
    if e:w=mul(w,powers[i][e],p)
   ms.append(w)
  num,den=linear(nv,ms),linear(dv,ms);assert den and num==mul(2,den,p)
  if uv is not None:
   aa=linear(proj[0],z);assert aa;inv=powmod(aa,(1<<d)-2,p)
   assert (mul(linear(proj[1],z),inv,p),mul(linear(proj[2],z),inv,p))==tuple(bits(c)for c in uv[k])
  zz.append(tuple(z));polys.append(p)
 assert len(set(zip(polys,zz)))==len(keys)
 if perm is not None:
  assert sorted(perm)==list(range(1,len(keys)+1))
  for k,(p,z)in enumerate(zip(polys,zz)):
   sz=tuple(norm([linear(row,z)for row in D],p));ix=perm[k]-1
   assert polys[ix]==p and zz[ix]==sz
  assert all(perm[perm[perm[perm[perm[k]-1]-1]-1]-1]==k+1 for k in range(len(keys)))
 return len(keys)
s=(H/'source_places_data.m').read_text();keys=lit(s,'SourcePlaceKeys');perm=lit(s,'SourceSigmaPermutation');orbs=lit(s,'SourceSigmaOrbits');indices=lit(s,'SourceJIndices')
assert check(keys,perm)==10 and len(orbs)==2 and all(len(o)==5 for o in orbs)
assert sorted(sum(orbs,[]))==list(range(1,11))
for orb in orbs:assert[perm[k-1]for k in orb]==orb[1:]+orb[:1]and len({indices[k-1]for k in orb})==1
assert[indices.count(k)for k in(1,2)]==[5,5]
targets={}
for i in(1,2,3):
 s=(H/f'target_places{i}_jbase_data.m').read_text();tk=lit(s,'TargetPlaceKeys');assert check(tk,lit(s,'TargetSigmaPermutation'))==5;assert all(k[0]==lit(s,'TargetJPolynomial')for k in tk);targets[str(i)]=5
result={'status':'PASS','source_places':10,'source_diamond_orbits':2,'source_orbit_lengths':[5,5],'canonical_quadrics_checked':276,'canonical_j_identity_checked_at_all_keys':True,'target_fibers_checked':targets}
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
print('NONIC71_FINITE_CANONICAL_KEYS_INDEPENDENT_PASS',json.dumps(result,sort_keys=True))
