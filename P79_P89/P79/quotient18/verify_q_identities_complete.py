#!/usr/bin/env python3
"""Read-only independent F2 q-expansion and linear-algebra replay for the genus27 model."""
from pathlib import Path
import re,ast,json,sys,hashlib
H=Path(__file__).resolve().parents[2]/'data/P79/quotient18';cs=(H/'canonical_data.m').read_text();js=(H/'jmap_q_data.m').read_text();N=600;mask=(1<<N)-1

def matrix(s,name):
 m=re.search(re.escape(name)+r':=Matrix\((?:GF\(2\)|Integers\(\)),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,name
 r,c=map(int,m.group(1,2));v=ast.literal_eval(m[3]);assert len(v)==r*c
 return r,c,[sum((v[i*c+j]&1)<<j for j in range(c))for i in range(r)]
def literal(name):return ast.literal_eval(re.search(re.escape(name)+r':=(.*?);',js,re.S)[1])
def rank(rows):
 b={}
 for v in rows:
  while v:
   k=v.bit_length()-1
   if k in b:v^=b[k]
   else:b[k]=v;break
 return len(b)
def mul(a,b):
 z=0
 while b:
  k=b&-b;z^=a<<(k.bit_length()-1);b^=k
 return z&mask

def act(v,M):
 z=0
 while v:
  k=v&-v;z^=M[k.bit_length()-1];v^=k
 return z
old=matrix(cs,'qcoeffs');q=matrix(js,'qcoeffs600');md=matrix(js,'Diamond3Matrix');rel=matrix(cs,'QuadRelations')
assert old[:2]==(18,130)and q[:2]==(18,600)and md[:2]==(18,18)and rel[:2]==(120,171)
assert [v&((1<<130)-1)for v in q[2]]==old[2]
assert rank([v&((1<<41)-1)for v in old[2]])==18
I=[1<<i for i in range(18)];p=I
for k in range(1,4):
 p=[act(v,md[2])for v in p]
 if k in[1,2]:assert p!=I
assert p==I
pairs=[tuple(map(int,x))for x in re.findall(r'<(\d+),\s*(\d+)>',re.search(r'pairs:=(.*?);',cs,re.S)[1])];assert len(pairs)==171 and set(pairs)=={(i,j)for i in range(1,19)for j in range(i,19)}
quad=[mul(q[2][i-1],q[2][j-1])for i,j in pairs];assert rank([v&((1<<130)-1)for v in quad])==51
for row in rel[2]:
 v=0
 while row:
  k=row&-row;v^=quad[k.bit_length()-1];row^=k
 assert v==0

def exps(name,d):
 out=literal(name)
 out=[v[0]if len(v)==1 and isinstance(v[0],list)else v for v in out]
 assert all(len(v)==18 and sum(v)==d and all(isinstance(a,int)and a>=0 for a in v)for v in out)
 assert len({tuple(v)for v in out})==len(out)
 return out
powers={}
for i,f in enumerate(q[2]):
 powers[i,0]=1
 for k in range(1,9):powers[i,k]=mul(powers[i,k-1],f)
def mon(v):
 r=1
 for i,k in enumerate(v):
  if k:r=mul(r,powers[i,k])
 return r
E2=list(map(mon,exps('Std2Exponents',2)));E8=list(map(mon,exps('Std8Exponents',8)))
assert len(E2)==51 and len(E8)==255 and rank(E2)==51 and rank(E8)==255 and rank(E2+E8)==291

def combine(vec,rows):
 assert len(vec)==len(rows)and set(vec)<={0,1};r=0
 for a,b in zip(vec,rows):
  if a:r^=b
 return r
G=combine(literal('Gvec'),E2);A=combine(literal('NumeratorVec'),E8);B=combine(literal('DenominatorVec'),E8)
assert G and A==G
Delta=2
for n in range(1,75):Delta=mul(Delta,1^(1<<(8*n))) if 8*n<N else Delta;Delta=mul(Delta,1^(1<<(16*n)))if 16*n<N else Delta
assert Delta==sum(1<<(m*m)for m in range(1,25,2)if m*m<N)
assert B==mul(Delta,G)and B
result={'N':79,'H_index':3,'genus':18,'q_basis_rank_through_Sturm40':18,'quadratic_relation_rank':120,'weight4_Sturm_bound':80,'quadratic_product_rank':51,'canonical_diamond3_exact_order':3,'degree8_canonical_product_rank':255,'auxiliary_G_space_dimension':15,'j_numerator_equals_E4cubed_times_G':True,'j_denominator_equals_Delta_times_G':True,'verified_q_precision':600,'weight16_Sturm_bound':320,'canonical_data_sha256':hashlib.sha256((H/'canonical_data.m').read_bytes()).hexdigest(),'jmap_data_sha256':hashlib.sha256((H/'jmap_q_data.m').read_bytes()).hexdigest(),'status':'Q_IDENTITIES_PASS; finite geometric j-map checks separately required'}
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
print('NONIC79_CANONICAL_Q_IDENTITIES_INDEPENDENT_PASS',json.dumps(result,sort_keys=True))
