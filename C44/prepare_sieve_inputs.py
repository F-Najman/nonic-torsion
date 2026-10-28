#!/usr/bin/env python3
from pathlib import Path
import argparse,json
ap=argparse.ArgumentParser();ap.add_argument('raw',type=Path);ap.add_argument('output',type=Path);a=ap.parse_args()
s=a.raw.read_text();dec=json.JSONDecoder();rows=[];i=0
while i<len(s):
 while i<len(s) and s[i].isspace():i+=1
 if i==len(s):break
 r,i=dec.raw_decode(s,i);rows.append(r)
assert rows[0]['type']=='header' and rows[-1]['type']=='complete'
n=rows[-1]['cases'];assert n==rows[0]['requested_cases'];assert rows[-1]['source_distinct_keys']==10*n
cases=[]
for k in range(1,n+1):
 j=[r for r in rows if r.get('case')==k and r['type']=='jet'];sides=[]
 for side in ['T','U']:
  rs=[r for r in j if r['side']==side];assert sum(r['degree']*r['multiplicity'] for r in rs)==54
  for r in rs:
   d=r['degree'];m=r['multiplicity'];assert len(r['field_polynomial'])==d+1 and len(r['canonical_jets'])==36
   assert all(len(z)==m and all(len(v)==d and set(v)<={0,1,2}for v in z)for z in r['canonical_jets'])
  sides.append(rs)
 cases.append((k,*sides))
def rec(r):return '[* '+str(r['degree'])+','+str(r['multiplicity'])+','+str(r['field_polynomial'])+','+str(r['canonical_jets'])+' *]'
with a.output.open('w')as f:
 f.write('JetCases:=[* *];\n')
 for k,ts,us in cases:f.write('Append(~JetCases,[* '+str(k)+', [* '+','.join(map(rec,ts))+' *], [* '+','.join(map(rec,us))+' *] *]);\n')
a.output.with_suffix('.normalized.jsonl').write_text(''.join(json.dumps(r,separators=(',',':'))+'\n' for r in rows))
print('C44_CANONICAL_SIEVE_INPUTS_PASS cases='+str(n))
