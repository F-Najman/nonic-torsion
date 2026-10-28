#!/usr/bin/env python3
from pathlib import Path
import json,time
D=Path(__file__).resolve().parent

def rank(rows):
 p={}
 for x in rows:
  while x:
   k=x.bit_length()-1
   if k in p:x^=p[k]
   else:p[k]=x;break
 return len(p)

def comps(parts,total,cap):
 if parts==1:
  if total<=cap:yield(total,)
  return
 for a in range(min(total,cap)+1):
  for b in comps(parts-1,total-a,cap):yield(a,)+b
s=(D/'data'/'FI_prime53_rows.txt').read_text().splitlines();p,ind,n,r,nrows,ncols=map(int,s[0].split())
assert(p,ind,n,r,nrows,ncols)==(53,13,2,2,57,1600)
assert len(s)==nrows+1 and all(len(x)==ncols and set(x)<={'0','1'}for x in s[1:])
rows=[int(x,2)for x in s[1:]];count=0;t=time.monotonic()
for d0 in range(1,10):
 for others in comps(ind-1,9-d0,min(d0,4)):
  sizes=(d0,)+others
  chosen=rows[:d0]
  for k in range(1,ind):chosen+=rows[9+(k-1)*4:9+(k-1)*4+sizes[k]]
  assert len(chosen)==9 and rank(chosen)==9,sizes
  count+=1
assert'NONIC53_ALLBAD_FI_MATRIX_EXPORT_COMPLETE'in(D/'logs'/'p53_allbad_FI_export.log').read_text()
out={'prime':53,'degree':9,'status':'ALLBAD2 BRANCH EXCLUDED ONLY','H_index':ind,'Tn':n,'q':3,'genus_quotient':40,'prefix_configurations':count,'rank_all':9,'binary_matrix_rank':rank(rows),'elapsed_seconds':time.monotonic()-t}
(D/'data'/'p53_allbad_verified.json').write_text(json.dumps(out,indent=2)+'\n')
print('NONIC53_ALLBAD_INDEPENDENT_BINARY_REPLAY_PASS',out)
