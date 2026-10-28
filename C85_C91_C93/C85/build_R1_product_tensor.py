#!/usr/bin/env python3
"""Exact R1xR1 -> R2 tensor via q1300, beyond weight4 Sturm1152."""
from pathlib import Path
import ast,re,json,sys,time,hashlib
H=Path(__file__).resolve().parents[1]/'data/C85';t0=time.monotonic();REPLAY='--readonly'in sys.argv;s=(H/'canonical_unit_stage1_data.m').read_text();m=re.search(r'qcoeffs:=Matrix\(GF\(2\),225,1300,(\[.*?\])\);',s,re.S);assert m;a=ast.literal_eval(m[1]);assert len(a)==225*1300;Q=[sum(int(a[1300*i+j])<<j for j in range(1300))for i in range(225)];mask=(1<<1300)-1

def mul(a,b):
 if a.bit_count()<b.bit_count():a,b=b,a
 z=0
 while b:
  k=b&-b;z^=a<<(k.bit_length()-1);b^=k
 return z&mask
piv={};basis=[];pairs=[];tensor=[[0]*225 for _ in range(225)]
for i in range(225):
 for j in range(i,225):
  original=mul(Q[i],Q[j]);x=original;c=0
  while x:
   k=x.bit_length()-1
   if k in piv:r,v=piv[k];x^=r;c^=v
   else:
    index=len(basis);piv[k]=(x,c^(1<<index));basis.append(original);pairs.append([i+1,j+1]);c=1<<index;break
  tensor[i][j]=tensor[j][i]=c
assert len(basis)==672 and len(piv)==672 and len(pairs)==672
# Independently reconstruct every product with the retained basis coefficients.
for i in range(225):
 for j in range(i,225):
  c=tensor[i][j];v=0
  while c:k=c&-c;v^=basis[k.bit_length()-1];c^=k
  assert v==mul(Q[i],Q[j])
out={'status':'PASS','dimensions':[225,672],'q_precision':1300,'Sturm_weight4':1152,'q_form_weight':2,'product_count':50625,'symmetric_products':25425,'basis_monomials':pairs,'basis_q':[hex(v)for v in basis],'product_1_1_into_2':[[hex(v)for v in row]for row in tensor],'native_data_sha256':hashlib.sha256(s.encode()).hexdigest()}
p=H/'canonical_ring_R1_tensor.json'
if REPLAY:assert json.loads(p.read_text())==out
else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
print('NONIC85_R1_PRODUCT_TENSOR_'+('READONLY_'if REPLAY else'')+'PASS dimR1=225 dimR2=672 products50625 q1300 Sturm1152 elapsed',time.monotonic()-t0,flush=True)
