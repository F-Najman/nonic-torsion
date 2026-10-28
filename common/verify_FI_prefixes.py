#!/usr/bin/env python3
"""Independent direct enumeration of every normalized cusp multiplicity tuple."""
from pathlib import Path
import sys,json
D=Path(sys.argv[1]).resolve();s=D.read_text().splitlines();vals=list(map(int,s[0].split()));p,ind,n,r,nrows,ncols=vals[:6]
assert len(s)==nrows+1 and nrows==9+4*(ind-1)
assert all(len(x)==ncols and set(x)<={'0','1'}for x in s[1:])
rows=[int(x,2)for x in s[1:]]
def rank(rr):
 pp={}
 for x in rr:
  while x:
   k=x.bit_length()-1
   if k in pp:x^=pp[k]
   else:pp[k]=x;break
 return len(pp)
def comps(k,total,cap):
 if k==1:
  if total<=cap:yield(total,)
  return
 for a in range(min(total,cap)+1):
  for bb in comps(k-1,total-a,cap):yield(a,)+bb
count=0
for d0 in range(1,10):
 for oo in comps(ind-1,9-d0,min(d0,4)):
  ss=(d0,)+oo;rr=rows[:d0]
  for k in range(1,ind):rr+=rows[9+4*(k-1):9+4*(k-1)+ss[k]]
  assert len(rr)==9 and rank(rr)==9,ss
  count+=1
res={'p':p,'degree':9,'H_index':ind,'Tn':n,'diamond_generator':r,'prefix_configurations':count,'rank_all':9,'full_binary_rank':rank(rows),'input_rows':nrows,'input_columns':ncols,'status':'PASS'}
if len(sys.argv)>2:Path(sys.argv[2]).write_text(json.dumps(res,indent=2)+'\n')
print('NONIC_INDEPENDENT_FI_PREFIX_PASS',json.dumps(res,sort_keys=True))
