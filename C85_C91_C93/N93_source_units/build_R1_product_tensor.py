#!/usr/bin/env python3
"""Exact R1xR1 -> R2 tensor via q1800, beyond weight4 Sturm1280."""
from pathlib import Path
import ast,re,json,sys,time,hashlib
H=Path(__file__).resolve().parents[1]/'data/N93_source_units';t0=time.monotonic();REPLAY='--readonly'in sys.argv;s=(H/'canonical_unit_data.m').read_text();m=re.search(r'qcoeffs:=Matrix\(GF\(2\),261,1800,(\[.*?\])\);',s,re.S);assert m;a=ast.literal_eval(m[1]);assert len(a)==261*1800;Q=[sum(int(a[1800*i+j])<<j for j in range(1800))for i in range(261)];mask=(1<<1800)-1

def mul(a,b):
 if a.bit_count()<b.bit_count():a,b=b,a
 z=0
 while b:
  k=b&-b;z^=a<<(k.bit_length()-1);b^=k
 return z&mask
piv={};basis=[];pairs=[];tensor=[[0]*261 for _ in range(261)]
for i in range(261):
 for j in range(i,261):
  original=mul(Q[i],Q[j]);x=original;c=0
  while x:
   k=x.bit_length()-1
   if k in piv:r,v=piv[k];x^=r;c^=v
   else:
    index=len(basis);piv[k]=(x,c^(1<<index));basis.append(original);pairs.append([i+1,j+1]);c=1<<index;break
  tensor[i][j]=tensor[j][i]=c
assert len(basis)==780 and len(piv)==780 and len(pairs)==780
# Independently reconstruct every product with the retained basis coefficients.
for i in range(261):
 for j in range(i,261):
  c=tensor[i][j];v=0
  while c:k=c&-c;v^=basis[k.bit_length()-1];c^=k
  assert v==mul(Q[i],Q[j])
out={'status':'PASS','dimensions':[261,780],'q_precision':1800,'Sturm_weight4':1280,'q_form_weight':2,'product_count':68121,'symmetric_products':34191,'basis_monomials':pairs,'basis_q':[hex(v)for v in basis],'product_1_1_into_2':[[hex(v)for v in row]for row in tensor],'native_data_sha256':hashlib.sha256(s.encode()).hexdigest()}
p=H/'canonical_ring_R1_tensor.json'
if REPLAY:assert json.loads(p.read_text())==out
else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
print('NONIC93_R1_PRODUCT_TENSOR_'+('READONLY_'if REPLAY else'')+'PASS dimR1=261 dimR2=780 products68121 q1800 Sturm1280 elapsed',time.monotonic()-t0,flush=True)
