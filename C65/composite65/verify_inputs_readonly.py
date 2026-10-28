#!/usr/bin/env python3
"""Bind completed C65 ordinary batches and exact winding/census provenance."""
from pathlib import Path
import json,hashlib,math
from unit_arithmetic import lit
H=Path(__file__).resolve().parents[1]/'data'/'composite65';R=H.parents[2];ss=[];tt=[];jj=[];fds=[]
for f,k in[(7,2),(8,2),(9,3)]:
 for i in range(1,k+1):
  P=H/f'hecke3_f{f}_source{i}';assert'NONIC65_ELLIPTIC_HECKE_JETS_STAGE2_COMPLETE'in(R/'C65'/'logs'/'composite65'/f'hecke3_f{f}_source{i}'/'run.log').read_text();s=(P/'elliptic_hecke3_data.m').read_text();ss+=lit(s,'SourceRecords');tt+=lit(s,'SourceHeckeRecords');jj+=lit(s,'SourceJPolynomials');fds.append(f)
s=(H/'elliptic_hecke_jets_stage2_data.m').read_text();assert ss==lit(s,'SourceRecords')and tt==lit(s,'SourceHeckeRecords')and jj==lit(s,'SourceJPolynomials')and fds==lit(s,'SourceDegrees')
s=(R/'C65'/'data'/'character_probes'/'N65_winding_data.m').read_text();HH={u for u in range(65)if math.gcd(u,65)==1}
for d,lev,w,o,k in lit(s,'ConstituentWindingFlags'):
 if not w:HH.intersection_update(k)
assert HH==set(lit(s,'FiniteDiamondUnits'))and len(HH)==48
print('NONIC65_ORDINARY_INPUT_BINDING_PASS completed_native_exits13 exact_combined_records=true finite_fullG24=true')
