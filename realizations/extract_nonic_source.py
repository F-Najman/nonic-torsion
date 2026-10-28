#!/usr/bin/env python3
"""Derive the exact nonic entries from the archived primary source."""
import hashlib
import json
from pathlib import Path
import re
import sympy as s

HERE=Path(__file__).resolve().parents[1]/'common'/'data'
source=HERE/'inputs'/'LowDegreePlaces'
x,y,t=s.symbols('x y t')
entries=[]
for lineno,line in enumerate(source.read_text().splitlines(),1):
    m=re.match(r'N = (\d+), degv = 9,.*?\[(.*) = 0, (.*) = 0\]',line)
    if not m:
        continue
    n=int(m[1]);f=s.sympify(m[2].replace('^','**'));g=s.sympify(m[3].replace('^','**'))
    assert s.degree(f,x)==9
    roots=s.solve(g,y);assert len(roots)==1
    xp=t;yp=roots[0].subs(x,t);ff=f.subs(x,t)
    def coeffs(p):
        return [str(c) for c in reversed(s.Poly(p,t).all_coeffs())]
    entries.append({'index':len(entries)+1,'N':n,'source_line':lineno,
                    'source_text':line,'f':coeffs(ff),'x':coeffs(xp),'y':coeffs(yp)})
assert len(entries)==19
out=HERE/'inputs'
data={'source':'https://www.math.fsu.edu/~hoeij/files/X1N/LowDegreePlaces',
      'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'examples':entries}
(out/'vanhoeij_nonic_input.json').write_text(json.dumps(data,indent=2)+'\n')
rows=['<'+str(e['N'])+','+', '.join('['+','.join(e[k])+']' for k in ('f','x','y'))+'>' for e in entries]
(out/'vanhoeij_nonic_input.m').write_text('NonicExamples:=[*\n'+',\n'.join(rows)+'\n*];\n')
print('NONIC_PRIMARY_SOURCE_EXTRACT_DONE examples=19 complete=true')
