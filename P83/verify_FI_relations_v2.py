#!/usr/bin/env python3
"""Read-only exact relation-kernel proof of all degree9 cusp restrictions."""
from pathlib import Path
import json,sys
P=Path(sys.argv[1]).resolve();s=P.read_text().splitlines();p,ind,n,r,nrows,ncols=map(int,s[0].split())
assert nrows==9+4*(ind-1)and len(s)==nrows+1
assert all(len(x)==ncols and set(x)<={'0','1'}for x in s[1:])
rows=[int(x,2)for x in s[1:]];piv={};kernel=[]
for i,x in enumerate(rows):
 v=1<<i
 while x:
  j=x.bit_length()-1
  if j in piv:a,b=piv[j];x^=a;v^=b
  else:piv[j]=(x,v);break
 if not x:kernel.append(v)
assert len(kernel)==nrows-len(piv)and len(kernel)<=20
minimum=None;checked=0;v=0;old=0
for j in range(1,1<<len(kernel)):
 gray=j^(j>>1);diff=gray^old;k=diff.bit_length()-1;v^=kernel[k];old=gray
 bits=v;first=bits&511;cost=first.bit_length();bits>>=9
 for _ in range(ind-1):cost+=(bits&15).bit_length();bits>>=4
 assert bits==0 and cost>9,(j,cost,v)
 minimum=cost if minimum is None else min(minimum,cost);checked+=1
count=0
for first in range(1,10):
 cap=min(first,4);dp=[1]+[0]*9
 for _ in range(ind-1):dp=[sum(dp[t-j]for j in range(min(cap,t)+1))for t in range(10)]
 count+=dp[9-first]
res={'p':p,'degree':9,'H_index':ind,'Tn':n,'diamond_generator':r,'input_rows':nrows,'input_columns':ncols,'full_binary_rank':len(piv),'relation_kernel_dimension':len(kernel),'checked_nonzero_relations':checked,'minimum_prefix_mass':minimum,'prefix_configurations':count,'rank_all':9,'method':'exhaustive complete relation kernel; every relation needs prefix mass>9','status':'PASS'}
if len(sys.argv)>2:Path(sys.argv[2]).write_text(json.dumps(res,indent=2)+'\n')
print('NONIC_INDEPENDENT_FI_RELATION_KERNEL_PASS',json.dumps(res,sort_keys=True))
