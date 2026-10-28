#!/usr/bin/env python3
"""Exact binary Laurent-series alignment of one serialized differential basis."""
from pathlib import Path
import ast,re,json,hashlib,time
H=Path(__file__).resolve().parents[1]/'data'/'quotient19';t0=time.monotonic();S=H/'plane_differential_certified_data.m';s=S.read_text();assert 'NONIC73_PLANE_DIFFERENTIAL_DATA_READY'in s

def lit(name):return ast.literal_eval(re.search(re.escape(name)+r':=(.*?);',s,re.S)[1].replace('<','(').replace('>',')'))
def mat(text,name):
 m=re.search(re.escape(name)+r':=Matrix\((?:GF\(2\)|Integers\(\)),(\d+),(\d+),(\[.*?\])\);',text,re.S);assert m,name
 r,c=map(int,m.group(1,2));a=ast.literal_eval(m[3]);assert len(a)==r*c
 return r,c,[sum((a[i*c+j]&1)<<j for j in range(c))for i in range(r)]
def val(a):assert a;return(a&-a).bit_length()-1
def mul(a,b,n):
 if a.bit_count()<b.bit_count():a,b=b,a
 z=0
 while b:
  k=b&-b;z^=a<<(k.bit_length()-1);b^=k
 return z&((1<<n)-1)
def square(a,n):
 z=0
 while a:
  k=a&-a;j=2*(k.bit_length()-1)
  if j<n:z|=1<<j
  a^=k
 return z
def inverse(a,n):
 assert a&1;b=1;k=1
 while k<n:
  k=min(2*k,n);b=mul(a,square(b,k),k)
 assert mul(a,b,n)==1
 return b
g=lit("genus");afd=lit("AF_degree");r,c,proj=mat(s,"ProjectionQSeries");assert r==3;va,vb,vc=map(val,proj);assert vb>va and vc>va
P=c-max(va,vb,vc)+min(vb-va,vc-va) # conservative common precision
# Each ratio actually has absolute precision1875; retain a conservative1872.
def quotient_pos(a,b):
 aa,bb=val(a),val(b);rel=1900-max(aa,bb);off=aa-bb;assert off>=0
 return(mul(a>>aa,inverse(b>>bb,rel),rel)<<off)&((1<<P)-1)
u,v=quotient_pos(proj[1],proj[0]),quotient_pos(proj[2],proj[0]);assert val(u)>0 and val(v)>0
base,gen=(u,v)if lit('AFBaseCoordinate')==1 else(v,u)
coefs=lit('DifferentialCoefficientData');assert len(coefs)==g and all(len(w)==afd for w in coefs)
polys={tuple(p)for w in coefs for c in w for p in c};maxdegree=max(map(len,polys))-1
powers=[1]
for k in range(1,maxdegree+1):powers.append(mul(powers[-1],base,P))
def poly(p):
 z=0
 for i,a in enumerate(p):
  assert a in[0,1]
  if a:z^=powers[i]
 return z
cache={p:poly(p)for p in polys};gp=[1]
for k in range(1,afd):gp.append(mul(gp[-1],gen,P))
du=0
for k in range(1,P,2):
 if u>>k&1:du|=1<<(k-1)
assert du
rows=[];worst=0
for w in coefs:
 terms=[]
 for j,(num,den)in enumerate(w):
  nn,dd=cache[tuple(num)],cache[tuple(den)];assert dd
  if nn==0:continue
  vn,vd=val(nn),val(dd);off=vn-vd;need=130-off;assert need>0 and need<=P-max(vn,vd),(off,need,P,vn,vd)
  b=mul(nn>>vn,inverse(dd>>vd,need),need);b=mul(b,gp[j],need);terms.append((off,b));worst=min(worst,off)
 low=min([0]+[off for off,b in terms]);z=0
 for off,b in terms:z^=b<<(off-low)
 z=mul(z,du,130-low)
 if low<0:assert z&((1<<(-low))-1)==0;z>>=-low
 row=(z<<1)&((1<<130)-1);assert row&1==0;rows.append(row)
old=mat((H/'canonical_data.m').read_text(),'qcoeffs')[2]
piv={}
for i,x in enumerate(rows):
 v=1<<i
 while x:
  k=(x&-x).bit_length()-1
  if k in piv:a,b=piv[k];x^=a;v^=b
  else:piv[k]=(x,v);break
assert len(piv)==g
change=[]
for x in old:
 v=0
 while x:
  k=(x&-x).bit_length()-1;assert k in piv
  a,b=piv[k];x^=a;v^=b
 change.append(v)
for target,comb in zip(old,change):
 z=0
 for i in range(g):
  if comb>>i&1:z^=rows[i]
 assert z==target
out=H/'canonical_change_basis_certified.m';values=[(v>>j)&1 for v in change for j in range(g)];out.write_text('CanonicalChangeBasis:=Matrix(GF(2),'+str(g)+','+str(g)+','+str(values)+');\n')
record={'input_sha256':hashlib.sha256(S.read_bytes()).hexdigest(),'genus':g,'plane_AF_degree':afd,'basis_dimension':g,'common_input_q_precision':P,'verified_q_precision':130,'differential_rank':g,'exact_match_to_saved_canonical_q_basis':True,'minimum_term_Laurent_offset':worst,'maximum_coefficient_polynomial_degree':maxdegree,'elapsed_seconds':time.monotonic()-t0,'status':'PASS; native exact field identities still required'}
(H/'canonical_change_basis_certified_verified.json').write_text(json.dumps(record,indent=2)+'\n');print('NONIC73_PLANE_DIFFERENTIAL_Q_ALIGNMENT_PASS',json.dumps(record,sort_keys=True))
