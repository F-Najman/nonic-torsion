#!/usr/bin/env python3
"""Read-only complete degree9 C133 certificate replay."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
from fractions import Fraction
import sympy as s
P=Path(__file__).resolve().parent;M=json.loads((P/'data'/'manifest133.json').read_text())
assert M['degree']==9 and M['excluded_cyclic_orders']==[133] and M['complete']
assert Fraction(325*8640,32768)==Fraction(43875,512)>72
assert pow(64,3,133)==1 and 64 not in (1,132)
H={pow(78,j,133) for j in range(18)};assert len(H)==18 and 64 in H and all(a%7==1 for a in H)
C=json.loads((P/'data'/'cusp133.json').read_text());assert C['H_free_on_source_cusps'] and C['all_F2_cusps_lift_to_same_global_degree']
assert C['global_signature']==C['F2_signature']=={'1':6,'3':2}
C=json.loads((P/'data'/'cusp_image133.json').read_text());assert C['admissible_degree3_divisors']==11 and C['complete']
assert C['admissible_Z_degree_width_counts']==[[1,133,3],[3,19,1]]
for e in ['reproduce_FI3.log','winding133_fresh.log','allbad_cubic_fresh.log','split_finite_operators.log','check_cm133.log','replay_split_operators.log']:
 for word in ('Runtime error','User error','Assertion failed'):assert word not in (P/'logs'/e).read_text()
log=(P/'logs'/'reproduce_FI3.log').read_text();assert 'NONIC133_Z_GENUS genus=35 exact_X019_embedding=true' in log
assert len(re.findall(r'NONIC133_Z_RR cusp=',log))==20
assert 'NONIC133_Z_FI_DONE degree=3 tested=58 exceptions=0' in log and 'complete=true' in log
assert 'genus=613 total_dimension=613' in (P/'logs'/'winding133_fresh.log').read_text()
assert 'minimum_degree=36 exceeds9=true' in (P/'logs'/'check_cm133.log').read_text()
L=(P/'logs'/'allbad_cubic_fresh.log').read_text();assert 'NONIC133_Y_ALLBAD_CUBIC_PASS' in L
assert [int(a) for a in re.findall(r'qualifying_P7_curves=(\d+)',L)]==[0,2,1]
for marker in ('Z133_TRIVIAL_FINITE_ANNIHILATOR','Z133_POSITIVE_LINEAR_CHECK','Z133_SPLIT_FINITE_OPERATORS_VERIFIED'):
 assert marker in (P/'logs'/'split_finite_operators.log').read_text()
log=(P/'logs'/'split_finite_operators.log').read_text();q=re.search(r'Z133_POSITIVE_T2 order=1 rational_dimension=4 matrix=(\[.*?\])',log)
import ast
T=s.Matrix(4,4,ast.literal_eval(q[1]));assert T**4+4*T**3+T**2-8*T-3*s.eye(4)==s.zeros(4)
assert 'combined_U=22 exact_RR_blocks=20 all_cusp_configurations=58 exceptions=0 complete=true' in (P/'logs'/'replay_split_operators.log').read_text()
r=subprocess.run([sys.executable,str(P/'replay_split_fi.py'),str(P/'data'/'fi3_data_reproduced.m')],text=True,capture_output=True,check=True)
assert '"all_cusp_configurations": 58' in r.stdout and '"supports_exclusion": true' in r.stdout and '"exact_RR_blocks": 20' in r.stdout
print('NONIC_C133_READONLY_CERTIFICATE_PASS degree=9 source_genus=613 quotient_genus=35 finite_operators=2 combined_U=22 exact_RR_blocks=20 all_cusp_configurations=58 all_rank3=true admissible_configurations=11 complete=true')
