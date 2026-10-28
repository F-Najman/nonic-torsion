#!/usr/bin/env python3
"""Read-only degree9 C81 proof-packet checks and binary matrix replay."""
from fractions import Fraction
from pathlib import Path
import hashlib,json,re,subprocess,sys
P=Path(__file__).resolve().parent
M=json.loads((P/'data'/'manifest81.json').read_text())
assert M['degree']==9 and M['excluded_cyclic_orders']==[81] and M['complete']
log=(P/'logs'/'fresh_native.log').read_text()
assert log.rstrip().endswith('NONIC81_NATIVE_DONE')
assert not re.search(r'Runtime error|User error|Assertion failed',log)
rows=[tuple(map(int,x)) for x in re.findall(r'NONIC81_WINDING order=(\d+) conductor=(\d+) dimension=(\d+) winding=(\d+) coefficient_degree=(\d+)',log)]
assert {(o,d,w,c) for o,_,d,w,c in rows}=={(1,4,4,1),(27,8,8,18),(9,6,6,6),(3,3,3,2)}
assert sum(d*c for o,f,d,w,c in rows)==190
assert re.findall(r'NONIC81_RATIONAL_TWO q=(\d+) integral_homology_rank=(\d+) joint_kernel_dimension=(\d+)',log)==[('5','380','8'),('7','380','0')]
R=subprocess.run([sys.executable,str(P/'degree8/replay_two_torsion_81.py'),str(P/'data/degree8/two_torsion_81_v2_matrices.m')],text=True,capture_output=True,check=True)
assert "'kernel_dimension': 0" in R.stdout and "'integral_homology_rank': 380" in R.stdout
index=Fraction(81*54,2)*Fraction(4,3);assert index==2916
cusps=27+2*9+6*3+18+27;assert 1+index/12-cusps/2==190
gon=Fraction(325,32768)*index;assert gon==Fraction(236925,8192)>3*9
assert M['principal_function_degree_bound']==27
assert M['requires_rational_cuspidal_torsion_equality'] is False
print('NONIC_C81_READONLY_CERTIFICATE_PASS degree=9 genus=190 rankQ=0 rational2torsion=0 Hecke_positive_degree=27 gonality=236925/8192 complete=true')
