#!/usr/bin/env python3
"""Exact Hasse and Waterhouse trace exclusion of order77 through F_(2^9)."""
import json,math
from pathlib import Path
H=Path(__file__).resolve().parent;rows=[]
for e in range(1,10):
 q=2**e;b=math.isqrt(4*q)
 multiples=[m for m in range(77,q+b+2,77) if abs(m-(q+1))<=b]
 if e%2:supersingular={0,2**((e+1)//2),-2**((e+1)//2)}
 else:supersingular={0,2**(e//2),-2**(e//2),2**(e//2+1),-2**(e//2+1)}
 traces=[q+1-m for m in multiples]
 allowed=[t for t in traces if t%2 or t in supersingular]
 assert not allowed
 rows.append({'e':e,'q':q,'Hasse_bound_floor':b,'multiples77_in_Hasse_interval':multiples,'candidate_traces':traces,'supersingular_traces':sorted(supersingular),'admissible_traces':allowed})
assert rows[-1]['multiples77_in_Hasse_interval']==[539] and rows[-1]['candidate_traces']==[-26]
(H/'data'/'finite_trace_certificate.json').write_text(json.dumps({'degree':9,'N':77,'complete':True,'source':'Waterhouse classification as stated in Derickx–Najman, arXiv2412.16016v2 Theorem4.6','rows':rows},indent=2)+'\n')
print('NONIC_77_FINITE_TRACE_PASS e=1..9 no_admissible_traces=true candidate_e9=-26 complete=true')
