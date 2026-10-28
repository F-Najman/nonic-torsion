#!/usr/bin/env python3
"""Independent small-target character and5961-support replay."""
from pathlib import Path
from itertools import combinations_with_replacement,product
from collections import Counter
import gzip,json,math,sys
import sympy as s
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parents[2]/'C35/data/presentation_v2'

def main():
 W=json.loads(gzip.decompress((P/'proven_pic_presentation.json.gz').read_bytes()));D=json.loads((P/'dual_characters.json').read_text());C=json.loads((P/'dual_classification.json').read_text())
 assert W['native_proven'] and W['factor_basis_bound']==9 and W['checkpoint_complete'] and D['complete'] and C['complete']
 assert W['known_Pic_order']==14529850982400 and math.prod(W['candidate_invariants'][:-1])==W['known_Pic_order']
 m=3592;assert len(W['factor_basis'])==m
 deg=[a['degree'] for a in W['factor_basis']];assert max(deg)==9
 # Newton identities, independently deriving every closed-place count through9.
 pc=W['frobenius_coefficients'];a=list(reversed(pc));assert a[0]==1 and len(pc)==51 and sum(pc)==W['known_Pic_order']
 power=[0];closed=[]
 for k in range(1,10):
  power.append(-k*a[k]-sum(a[j]*power[k-j] for j in range(1,k)))
  pts=3**k+1-power[k];mub=pts-sum(d*closed[d-1] for d in range(1,k) if k%d==0);assert mub%k==0;closed.append(mub//k)
 assert closed==[deg.count(i) for i in range(1,10)]==[12,0,12,3,36,128,348,813,2240]
 B2=D['character_matrix16'];B5=D['character_matrix5'];r=len(B2);assert r==D['character_generators16'] and len(B5)==2 and all(len(v)==m for v in B2+B5)
 base=D['base_FB_index']-1;assert deg[base]==1 and C['base_FB_index']==base+1 and all(v[base]==0 for v in B2+B5)
 for row in W['relations']:
  assert len(row)==m;sparse=[(i,x) for i,x in enumerate(row) if x];assert sum(x*deg[i] for i,x in sparse)==0
  assert all(sum(x*b[i] for i,x in sparse)%16==0 for b in B2)
  assert all(sum(x*b[i] for i,x in sparse)%5==0 for b in B5)
 B2M=s.Matrix(B2);B5M=s.Matrix(B5)
 def lattice(B,mod,inds=None):
  V=B if inds is None else B[:,[i-1 for i in inds]]
  H=hermite_normal_form(V.row_join(mod*s.eye(B.rows)));return H
 I2=lattice(B2M,16);I5=lattice(B5M,5)
 assert 16**r//int(abs(I2.det()))==16384 and 5**2//int(abs(I5.det()))==25
 cusps=C['cusp_FB_indices'];assert len(cusps)==17 and len(set(cusps))==17 and sorted(deg[i-1] for i in cusps)==[1]*12+[4]*3+[6]*2
 H2=lattice(B2M,16,cusps);H5=lattice(B5M,5,cusps);assert 16**r//int(abs(H2.det()))==8 and 25//int(abs(H5.det()))==5
 assert (16384*25)//(8*5)==10240 and W['known_Pic_order']//10240==13*109148520
 H2i,H5i=H2.inv(),H5.inv()
 A2=(16*H2i).tolist();A5=(5*H5i).tolist();assert all(x.q==1 for row in A2+A5 for x in row);A2=[[int(x) for x in row] for row in A2];A5=[[int(x) for x in row] for row in A5]
 def in_small(A,v,mod):return all(sum(a*b for a,b in zip(row,v))%mod==0 for row in A)
 def member(inv,v):return all(x.q==1 for x in inv*s.Matrix(v))
 nc=C['noncuspidal_FB_indices'];assert len(nc)==9 and list(map(len,nc))==[0,0,12,0,36,126,348,813,2240]
 assert sorted(cusps+[i for row in nc for i in row])==list(range(1,m+1))
 assert all(deg[i-1]==d for d,row in enumerate(nc,1) for i in row)
 known=set(C['known_octic_FB_indices']);assert len(known)==6 and known<=set(nc[7])
 tests={};unexpected=[];gotknown=set()
 for T in C['tests']:
  key=(tuple(T['degrees']),tuple(T['indices']));assert key not in tests;tests[key]=T
 expected=[]
 for d in (3,5,6,7,8,9):expected += [((d,),(i,)) for i in range(1,len(nc[d-1])+1)]
 expected += [((3,3),ij) for ij in combinations_with_replacement(range(1,13),2)]
 expected += [((d,3),(i,j)) for d in (5,6) for i in range(1,len(nc[d-1])+1) for j in range(1,13)]
 expected += [((3,3,3),ijk) for ijk in combinations_with_replacement(range(1,13),3)]
 assert len(expected)==5961 and set(expected)==set(tests)
 for typ,inds in expected:
  T=tests[typ,inds];ids=[nc[d-1][j-1]-1 for d,j in zip(typ,inds)]
  v2=[sum(row[i] for i in ids)%16 for row in B2];v5=[sum(row[i] for i in ids)%5 for row in B5]
  assert T['image16']==v2 and T['image5']==v5
  inh=in_small(A2,v2,16) and in_small(A5,v5,5);isknown=typ==(8,) and ids[0]+1 in known
  assert T['in_H']==inh and T['known_octic_place']==isknown
  if inh:
   if isknown:gotknown.add(ids[0]+1)
   else:unexpected.append({'degrees':typ,'indices':inds})
 report={'level':35,'native_proven_Pic_order':W['known_Pic_order'],'factor_basis_places':3592,'generation_bound':9,'primary_image_orders':[16384,25],'rational_primary_image_orders':[8,5],'rational_image_index':10240,'all_relation_rows_checked':len(W['relations']),'supports':5961,'known_octic_places':len(gotknown),'unexpected':unexpected,'supports_exclusion':not unexpected and gotknown==known,'complete':True}
 print('NONIC35_PRIMARY_INTEGER_REPLAY '+json.dumps(report))
 if '--write-report' in sys.argv:(P/'dual_integer_replay.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
