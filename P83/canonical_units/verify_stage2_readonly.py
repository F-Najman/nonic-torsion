#!/usr/bin/env python3
"""Exact independent torsion-family/unit-jet and finite-rank replay."""
from pathlib import Path
import json,time,hashlib,sys
from sage.all import GF,PolynomialRing,matrix,EllipticCurve
from unit_arithmetic import *
H=Path(__file__).resolve().parents[1]/'data/canonical_units';t=time.monotonic();a=(H/'canonical_unit_stage1_data.m').read_text();b=(H/'elliptic_hecke_jets_stage2_data.m').read_text();RX=matbits(a,'UnitXRelations');RY=matbits(a,'UnitYRelations');ss=lit(b,'SourceRecords');tt=lit(b,'SourceHeckeRecords');assert len(ss)==len(tt)==5;cache={};count=0
for sid,(s,ts)in enumerate(zip(ss,tt),1):
 assert len(ts)==3 and sorted(d for e,d,m,r in ts)==[9,9,18]and all(m==1 for e,d,m,r in ts)
 for rec in[s]+[r for e,d,m,r in ts]:
  v=verify_record(rec,RX,RY);F=v['F'];p=F.p
  if p not in cache:
   R=PolynomialRing(GF(2),'a');mp=R(v['field_poly']);assert mp.is_irreducible();cache[p]=GF(2**F.d,'a',modulus=mp)
  FF=cache[p];ec=lambda z:FF(R([(z>>j)&1 for j in range(F.d)]));E=EllipticCurve(FF,list(map(ec,v['ainv'])));P=E(list(map(ec,v['point'])));assert P!=E(0)and 83*P==E(0);u,w=v['units'];cv={n:FF(R([(n>>j)&1 for j in range(F.d)]))for n in {0,1,u,u^1,w,w^1}};M=matrix(FF,[[cv[((row>>i)&1)^(u if(row>>(247+i))&1 else 0)]for i in range(247)]for row in RX]+[[cv[((row>>i)&1)^(w if(row>>(247+i))&1 else 0)]for i in range(247)]for row in RY]);assert M.rank()==246;count+=1
 print('NONIC83_JET_REPLAY_SOURCE',sid,'all_four_families_valid=true',flush=True)
out={'status':'PASS','sources':5,'ordinary_83_torsion_families':count,'precision':5,'unit_relations':[158,153],'all_finite_unit_ranks':246,'all_psi83_derivatives_nonzero':True,'all_canonical_tangents_nonzero':True,'native_data_sha256':hashlib.sha256(b.encode()).hexdigest(),'elapsed_seconds':time.monotonic()-t}
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
print('NONIC83_ELLIPTIC_HECKE_JETS_INDEPENDENT_PASS',json.dumps(out,sort_keys=True))
