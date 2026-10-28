#!/usr/bin/env python3
"""Select small exact good reductions for an independent direct point count."""
import json,re
from math import gcd
from pathlib import Path
H=Path(__file__).resolve().parents[1]
E=json.loads((H/'common/data/inputs/vanhoeij_nonic_input.json').read_text())['examples']
log=(H/'realizations/logs/verify_nonic_examples.log').read_text()
allrec=re.findall(r'NONIC_GOOD_REDUCTION index=(\d+) p=(\d+) residue_degree=(\d+) order=(\d+) gcd=(\d+)',log)
records=[]
for e in E:
    opts=sorted(set(tuple(map(int,z[1:4])) for z in allrec if int(z[0])==e['index']),key=lambda x:(x[0]**x[1],x))
    bound=0;selected=[]
    for p,d,n in opts:
        g=gcd(bound,n)
        if bound==g:continue
        selected.append([p,d,n]);bound=g
        if bound==e['N']:break
    assert bound==e['N']
    # Remove a reduction when its deletion preserves the exact bound.
    for i in range(len(selected)-1,-1,-1):
        b=0
        for j,rec in enumerate(selected):
            if i!=j:b=gcd(b,rec[2])
        if b==e['N']:selected.pop(i)
    records.append({'index':e['index'],'N':e['N'],'reductions':selected})
(H/'realizations/data/selected_reductions.json').write_text(json.dumps(records,indent=2)+'\n')
rows=['<'+str(e['index'])+','+'['+','.join('<'+','.join(map(str,r))+'>' for r in e['reductions'])+']>' for e in records]
(H/'realizations/data/selected_reductions.m').write_text('SelectedReductions:=[*\n'+',\n'.join(rows)+'\n*];\n')
for r in records:print(r)
print('NONIC_SMALL_REDUCTIONS_SELECTED count=19 complete=true')
