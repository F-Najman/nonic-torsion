#!/usr/bin/env python3
"""Exact binary Laurent-series alignment of one serialized differential basis."""
from pathlib import Path
import ast,re,json,hashlib,time
H=Path(__file__).resolve().parents[1]/'data'/'quotient45';t0=time.monotonic();S=H/'plane_differential_certified_data.m';s=S.read_text();assert 'NONIC67H3_PLANE_DIFFERENTIAL_DATA_READY'in s

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
 aa,bb=val(a),val(b);rel=c-max(aa,bb);off=aa-bb;assert off>=0
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
  vn,vd=val(nn),val(dd);off=vn-vd;need=2200-off;assert need>0 and need<=P-max(vn,vd),(off,need,P,vn,vd)
  b=mul(nn>>vn,inverse(dd>>vd,need),need);b=mul(b,gp[j],need);terms.append((off,b));worst=min(worst,off)
 low=min([0]+[off for off,b in terms]);z=0
 for off,b in terms:z^=b<<(off-low)
 z=mul(z,du,2200-low)
 if low<0:assert z&((1<<(-low))-1)==0;z>>=-low
 row=(z<<1)&((1<<2200)-1);assert row&1==0;rows.append(row)
import sys
REPLAY='--readonly'in sys.argv;WORK=2200
change=mat((H/'canonical_change_basis_certified.m').read_text(),'CanonicalChangeBasis')[2]
def xor_rows(v,rs):
 out=0
 while v:
  b=v&-v;out^=rs[b.bit_length()-1];v^=b
 return out
qext=[xor_rows(v,rows)for v in change];js=(H/'jmap_q_data.m').read_text();saved=mat(js,'qcoeffs1400')[2];assert len(qext)==45 and[a&((1<<1400)-1)for a in qext]==saved
D=mat(js,'Diamond2Matrix')[2];qsig=[xor_rows(v,qext)for v in D]
def jlit(n):return ast.literal_eval(re.search(re.escape(n)+r':=(.*?);',js,re.S)[1].replace('<','(').replace('>',')'))
exps=jlit('Std10Exponents');nv=jlit('NumeratorVec');dv=jlit('DenominatorVec');assert len(exps)==len(nv)==len(dv)==836

def evaluate_j_pair(qs):
 ps=[[1]for _ in qs]
 for i,a in enumerate(qs):
  for _ in range(10):ps[i].append(mul(ps[i][-1],a,WORK))
 num=den=0
 for e,a,b in zip(exps,nv,dv):
  if not a and not b:continue
  v=1
  for i,k in enumerate(e):
   if k:v=mul(v,ps[i][k],WORK)
  if a:num^=v
  if b:den^=v
 return num,den
num,den=evaluate_j_pair(qext);sn,sd=evaluate_j_pair(qsig);assert num and den and sn and sd
assert mul(sn,den,WORK)==mul(num,sd,WORK)
# The homogeneous cross-product is a section of20K. Its q-series uses
# (dq/q)^20, so zero throughq2200 gives vanishingorder at least2180.
assert WORK-20>20*(2*45-2)
out={'status':'PASS','q_precision':WORK,'genus':45,'native_q_prefix_verified':1400,'canonical_section_degree':20,'canonical_section_line_bundle_degree':1760,'proved_vanishing_order_at_least':2180,'sigma_preserves_j':True,'canonical_q_series':[hex(a)for a in qext],'plane_input_sha256':hashlib.sha256(S.read_bytes()).hexdigest(),'change_basis_sha256':hashlib.sha256((H/'canonical_change_basis_certified.m').read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t0}
f=H/'canonical_extended_q_sigma_verified.json'
if REPLAY:
 old=json.loads(f.read_text())
 for k in out:
  if k!='elapsed_seconds':assert old[k]==out[k],k
else:f.write_text(json.dumps(out,indent=2)+'\n')
print('NONIC67H3_EXTENDED_Q_SIGMA_CERTIFICATE_PASS precision2200 sigma_j_invariant=true section20K_degree1760 vanishing>=2180 elapsed',time.monotonic()-t0)
