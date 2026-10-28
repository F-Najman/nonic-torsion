#!/usr/bin/env python3
"""Bind the independently reconstructed C91 data to completed original producers."""
from pathlib import Path
import json,hashlib
from unit_arithmetic import lit
H=Path(__file__).resolve().parents[1]/'data/C91';R=H.parents[2];U=R/'C85_C91_C93/data/N91_source_units';UL=R/'C85_C91_C93/logs/N91_source_units'
assert(H/'canonical_unit_stage1_data.m').read_bytes()==(U/'canonical_unit_data.m').read_bytes()
a=(H/'elliptic_hecke_jets_stage2_data.m').read_text();ss=[];tt=[];jj=[]
for i in range(1,5):
 P=U/f'hecke3_source{i}';assert'NONIC91_ELLIPTIC_HECKE_JETS_STAGE2_COMPLETE'in(UL/f'hecke3_source{i}'/'run.log').read_text();b=(P/'elliptic_hecke3_data.m').read_text();ss+=lit(b,'SourceRecords');tt+=lit(b,'SourceHeckeRecords');jj+=lit(b,'SourceJPolynomials')
assert ss==lit(a,'SourceRecords')and tt==lit(a,'SourceHeckeRecords')and jj==lit(a,'SourceJPolynomials')and len(ss)==4
s=(R/'C85_C91_C93/data/character_probes/N91_winding_data.m').read_text();flags=lit(s,'ConstituentWindingFlags');HH={u for u in range(91)if __import__('math').gcd(u,91)==1}
for d,lev,w,o,k in flags:
 if not w:HH.intersection_update(k)
assert HH==set(lit(s,'FiniteDiamondUnits'))and 3 in HH and 27 in HH and len(HH)==12
print('NONIC91_INPUT_BINDING_PASS completed_native_exits10 exact_combined_records=true finite_diamond3_and27=true')
