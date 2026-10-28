#!/usr/bin/env python3
"""Read-only complete coverage audit, including explicit degree-nine transfer conditions."""
from pathlib import Path
import json,re,hashlib,sys
H=Path(__file__).resolve().parent/'logs'
# --rerun: audit Z/46 on the complete rerun of the whole degree-9 sieve instead of the original combination of runs
RERUN='--rerun'in sys.argv

def parts(n,b=None):
 if n==0:yield ();return
 for a in range(min(n,b or n),0,-1):
  for rest in parts(n-a,a):yield(a,)+rest
assert sum(sum(1 for _ in parts(i))for i in range(1,10))==96
ind=(H.parents[1]/'common/logs/independent_counts_46_50_51.log').read_text()
for n,old,new,files in [(46,{(4,):5,(5,):4,(6,):5,(7,):24,(8,):71,(4,4):275},{(9,):178,(5,4):220},['N46_p3_single_completed.log','N46_p3_mixed_quotient.log']),(50,{(6,):6,(7,):20,(8,):107},{(9,):204},['N50_p3_quotient.log'])]:
 if RERUN and n==46:old,new,files={},{**old,**new},[]
 texts=[(H/'degree8'/fn).read_text()for fn in files];s=(H/('N46_p3_quotient_D9_rerun.log'if RERUN and n==46 else f'N{n}_p3_good9.log')).read_text()
 assert f'SUMMARY N={n} M=1 p=3 D=9 shard=0/1 total_survivors=0'in s and s.rstrip().endswith('SIEVE_DONE'if RERUN and n==46 else'SIEVE_RESUME_DONE FromDegree=9 MaxPart=9')
 for fn,t in zip(files,texts):
  if 'single_completed'in fn:assert t.rstrip().endswith('SINGLE_BRANCH_CHECKPOINT_DONE level=46 p=3 D=8 completed_degrees=4,5,6,7,8')
  else:assert f'SUMMARY N={n} M=1 p=3 D=8 shard=0/1 total_survivors=0'in t and t.rstrip().endswith(f'SIEVE_RESUME_DONE FromDegree={8 if n==46 else 6} MaxPart={4 if n==46 else 8}')
 seen={}
 for t,expected in [( '\n'.join(texts),old),(s,new)]:
  assert not re.search(r'Runtime error|User error|System error|Assertion failed|syntax error|failed memory request|^ERROR',t,re.I|re.M)
  got={}
  for m in re.finditer(r'RESULT type=\[([^\]]+)\] total=(\d+) survivors=(\d+)',t):
   part=tuple(map(int,m[1].replace(',',' ').split()));assert part not in got and m[3]=='0';got[part]=int(m[2])
   indices=[int(x[1])for x in re.finditer(r'SMART_CASE_RESULT type=\[\s*'+r'\s*,\s*'.join(map(str,part))+r'\s*\] index=(\d+) killed=true',t)]
   if indices or expected is new or (n==46 and len(part)==1):assert sorted(indices)==list(range(1,int(m[2])+1)),(n,part)
   else:assert 'SIEVE_RESUME_DONE' in t and int(m[2])==expected[part]
  assert got==expected,(n,got,expected);seen.update(got)
 c=list(map(int,re.search(r'NONCUSP counts deg 1\.\.9: \[([^\]]+)\]',s)[1].replace(',',' ').split()))
 ci=list(map(int,re.search(rf'DEGREE9_FAST_COUNT m=1 n={n} p=3 .*?closed_place_counts=\[([^\]]+)\]',ind)[1].replace(',',' ').split()));assert c==ci
 possible={p for d in range(1,10)for p in parts(d)if all(c[a-1]for a in p)};assert possible==set(seen)
 # Every old class rejection is still valid: the larger target-cusp complements all have quotient class zero.
 assert 'complement_counts=[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]'in s and 'predicted F_p cusp degrees match model: true'in s
 maxcomp=9-min(c.index(x)+1 for x in c if x)
 expectedq=[5,7,11,13,17,19]if n==46 else[7,11,13,17,19]
 for k in range(1,maxcomp+1):
  q=list(map(int,re.search(rf'CC complement degree <= {k}:.*?usable q from list: \[([^\]]+)\]',s)[1].replace(',',' ').split()));assert q==expectedq
 print('NONIC_LOWER_AND_NEW_COVERAGE_PASS',json.dumps({'N':n,'degree':9,'counts':c,'support_types':[{'type':list(k),'tests':v}for k,v in sorted(seen.items())],'old_tests':sum(old.values()),'new_tests':sum(new.values()),'total_tests':sum(seen.values()),'max_cusp_complement':maxcomp,'all_96_partitions_accounted_for':True,'zero_survivors':True}))
print('NONIC46_50_READONLY_COVERAGE_PASS')
