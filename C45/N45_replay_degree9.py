#!/usr/bin/env python3
"""Check all degree9 cusp restrictions except pure9 rational cusps.
Uses immutable degree-independent jets through weighted degree8, without
mixing differential bases. Pure9cusps require a separate fresh native test.
"""
from pathlib import Path
import re,ast,json,time,hashlib,sys
H=Path(__file__).resolve().parent;P=H/'data/degree8/N45_completion_jets.m';s=P.read_text();t0=time.monotonic()
assert s.count('FI_JET_DATA_READY')==1
weights=ast.literal_eval(re.search(r'ds:=(.*?);',s)[1]);assert sorted(weights)==[1]*12+[2]*4+[4]*5+[6]*2
dim=int(re.search(r'udim:=(\d+)',s)[1]);assert dim==39
jets={}
for i,m,r,c,txt in re.findall(r'localjets\[<(\d+),(\d+)>\]:=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S):
 i,m,r,c=map(int,(i,m,r,c));v=ast.literal_eval(txt);assert r==39 and c==m*weights[i-1]and len(v)==r*c
 cols=tuple(sum(v[k*c+j]<<k for k in range(r))for j in range(c));jets[i-1,m]=cols
for i,w in enumerate(weights):
 for m in range(1,8//w+1):
  assert (i,m)in jets
  if m>1:assert jets[i,m][:w*(m-1)]==jets[i,m-1]
maxcols=[jets[i,8//w]for i,w in enumerate(weights)];size=len(weights);ways=[[0]*10 for _ in range(size+1)];ways[size][0]=1
for i in range(size-1,-1,-1):
 for r in range(10):ways[i][r]=sum(ways[i+1][r-m*weights[i]]for m in range(min(8//weights[i],r//weights[i])+1))
expected=ways[0][9]
basis=[0]*39;vec=[0]*size;count=0;exceptions=[];exception_count=0

def visit(i,remaining,rank):
 global count,exception_count
 if remaining==0:
  count+=1
  if rank<9:
   exception_count+=1
   if len(exceptions)<10000:exceptions.append({'rank':rank,'vector':vec.copy()})
  if count%250000==0:print('N45_DEGREE9_CUSP_PROGRESS',count,exception_count,round(time.monotonic()-t0,2),flush=True)
  return
 if i==size or ways[i][remaining]==0:return
 w=weights[i];maximum=min(8//w,remaining//w);added=[];current=rank
 for m in range(maximum+1):
  if m:
   for col in maxcols[i][(m-1)*w:m*w]:
    while col:
     p=col.bit_length()-1
     if basis[p]:col^=basis[p]
     else:basis[p]=col;added.append(p);current+=1;break
  vec[i]=m
  if remaining-m*w and 9-(remaining-m*w)<=8:assert current==9-(remaining-m*w)
  visit(i+1,remaining-m*w,current)
 vec[i]=0
 for p in added:basis[p]=0
visit(0,9,0)
assert count==expected
result={'degree':9,'source_order':45,'geometry_file_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'tested_divisors':count,'exceptions_count':exception_count,'exceptions':exceptions,'pure9_rational_cusps_separate_control':12,'all_cuspidal_divisors':count+12,'complete_enumeration':True,'elapsed_seconds':time.monotonic()-t0,'status':'FULL_RANK'if not exception_count else'EXCEPTIONS REQUIRE FURTHER ARGUMENT'}
if '--readonly' not in sys.argv:(H/'data'/'N45_degree9_FI_replay.json').write_text(json.dumps(result,indent=2)+'\n')
else:
 stored=json.loads((H/'data'/'N45_degree9_FI_replay.json').read_text());assert all(stored[k]==v for k,v in result.items()if k!='elapsed_seconds')
print('N45_DEGREE9_CUSP_REPLAY_DONE tested='+str(count)+' exceptions='+str(exception_count)+' complete=true')
