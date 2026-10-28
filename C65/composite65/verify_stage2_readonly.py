#!/usr/bin/env python3
"""All C65 ordinary torsion families and unit-recovered canonical jets."""
from pathlib import Path
import json,sys,time,hashlib,os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from sage.all import GF,PolynomialRing,EllipticCurve,matrix
from unit_arithmetic import *
H=Path(__file__).resolve().parents[1]/'data'/'composite65';t=time.monotonic();a=(H.parent/'source_units'/'canonical_unit_data.m').read_text();b=(H/'elliptic_hecke_jets_stage2_data.m').read_text();RX=matbits(a,'UnitXRelations');RY=matbits(a,'UnitYRelations');ss=lit(b,'SourceRecords');tt=lit(b,'SourceHeckeRecords');fds=lit(b,'SourceDegrees');assert fds==[7,7,8,8,9,9,9]and len(ss)==len(tt)==7;cache={};count=0;patterns=[]
for sid,(f,s,ts)in enumerate(zip(fds,ss,tt),1):
 assert all(m==1 and d==f*e for e,d,m,r in ts)and sum(d*m for e,d,m,r in ts)==4*f;patterns.append([d for e,d,m,r in ts])
 for rec in[s]+[r for e,d,m,r in ts]:
  v=verify_record(rec,RX,RY);F=v['F'];p=F.p
  if p not in cache:
   R=PolynomialRing(GF(2),'a');mp=R(v['field_poly']);assert mp.is_irreducible();cache[p]=GF(2**F.d,'a',modulus=mp)
  FF=cache[p];ec=lambda z:FF(R([(z>>j)&1 for j in range(F.d)]));E=EllipticCurve(FF,list(map(ec,v['ainv'])));P=E(list(map(ec,v['point'])));assert P!=E(0)and 65*P==E(0)and 5*P!=E(0)and 13*P!=E(0);u,w=v['units'];cv={n:FF(R([(n>>j)&1 for j in range(F.d)]))for n in{0,1,u,u^1,w,w^1}};M=matrix(FF,[[cv[((row>>i)&1)^(u if(row>>(121+i))&1 else 0)]for i in range(121)]for row in RX]+[[cv[((row>>i)&1)^(w if(row>>(121+i))&1 else 0)]for i in range(121)]for row in RY]);assert M.rank()==120;count+=1
 print('NONIC65_ORDINARY_JET_SOURCE_PASS',sid,'degree',f,flush=True)
assert patterns==[[28],[28],[16,16],[16,16],[36],[36],[36]]and count==16
out={'status':'PASS','source_degrees':fds,'ordinary65_torsion_families':16,'precision':5,'unit_relations':[69,68],'all_finite_unit_ranks':120,'all_psi65_derivatives_nonzero':True,'all_canonical_tangents_nonzero':True,'target_degree_patterns':patterns,'native_combined_data_sha256':hashlib.sha256(b.encode()).hexdigest()}
p=H/'stage2_independent_verified.json'
if '--readonly'in sys.argv:assert json.loads(p.read_text())==out
else:p.write_text(json.dumps(out,indent=2)+'\n')
print('NONIC65_ORDINARY_STAGE2_PASS elapsed',time.monotonic()-t)
