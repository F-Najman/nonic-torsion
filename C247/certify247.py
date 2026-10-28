#!/usr/bin/env python3
"""Read-only complete C247 exclusion, including independent norm replay."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,re,subprocess,sys
from math import gcd
import sympy as sp
P=Path(__file__).resolve().parent;M=json.loads((P/'data'/'manifest247.json').read_text())
assert M['degree']==9 and M['excluded_cyclic_orders']==[247] and M['complete']
log=(P/'logs'/'winding247.log').read_text()
assert log.rstrip().endswith('WINDING_CHARACTER_PROBE_DONE') and not re.search('Runtime error|User error|Assertion failed',log)
rows=[tuple(map(int,z)) for z in re.findall(r'WINDING N=247 order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+)',log)]
assert sum(d*f for o,c,d,w,f in rows)==2305
assert 'admissible=[ 12, 56, 68, 69, 75, 87, 88, 103, 144, 159, 160, 172, 178, 179, 191, 235 ]' in log
log=(P/'logs'/'positive_polynomial.log').read_text()
rows=[tuple(map(int,z)) for z in re.findall(r'Z247_WINDING order=(\d+) conductor=(\d+) dimension=(\d+) winding=(\d+) field_degree=(\d+)',log)]
assert {(o,d,w,f) for o,c,d,w,f in rows}=={(1,21,14,1),(2,20,20,1),(3,22,22,2),(6,22,22,2)}
assert sum(d*f for o,c,d,w,f in rows)==129
assert 129-21==108 and 2*22+2*22==88
n=247;phi=int(sp.totient(n));index=Fraction(n*phi,2)*Fraction(14,13)*Fraction(20,19)
assert index==30240
cusps=sum(int(sp.totient(d))*int(sp.totient(n//d)) for d in sp.divisors(n))//2
assert cusps==432 and 1+index/12-cusps/2==2305
assert Fraction(325,32768)*index==Fraction(307125,1024)>72
assert (2305-1)//18+1==129 and index/18==1680
assert Fraction(325,32768)*(index/18)==Fraction(34125,2048)>15
assert pow(144,3,n)==1 and all(pow(144,k,n) not in [1,n-1] for k in [1,2])
H={pow(40,k,n) for k in range(18)};assert len(H)==18 and all(h%13==1 for h in H) and 144 in H
assert pow(2,6,13)==12 and all(pow(2,k,13) not in [1,12] for k in range(1,6))
C=json.loads((P/'data'/'cusp247.json').read_text());assert C['H_free_on_source_cusps'] and C['source_cusp_count']==432 and C['quotient_cusp_count']==24
cm={int(a):int(b) for a,b,_ in re.findall(r'\[\*\s*(\d+),\s*(\d+),\s*(\d+),',(P.parent/'common/data/dcm_list_all_min_orders_X1_10k.m').read_text())};assert cm[247]==72>9
r=subprocess.run([sys.executable,str(P/'replay_norms.py')],text=True,capture_output=True,check=True)
assert r.stdout.count('INDEPENDENT_NORM_PASS')==3 and 'PARITY_REPLAY_COMPLETE' in r.stdout
print('NONIC_C247_READONLY_CERTIFICATE_PASS degree=9 source_genus=2305 quotient_genus=129 source_principal_degree=72 target_principal_degree=15 connected_A_dimension=108 rational_A_2torsion=0 three_exact_norms_replayed=true complete=true')
