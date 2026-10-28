#!/usr/bin/env python3
"""Independent coverage audit for the completed nonic C38 sieve branches."""
from pathlib import Path
from collections import Counter
from math import comb,prod
import ast,json,re,sys
H=Path(__file__).resolve().parent/'logs'
names=sys.argv[1:] or ['N38_p3_quotient_D9_rerun.log'];logs=[(H/n).read_text()for n in names]
nc=[0,0,3,0,72,87,216,801,2127];reps=[0,0,1,0,8,11,24,89,239]
for text in logs:
 assert not re.search(r'Runtime error|Assertion failed|User error|ERROR|SURVIVOR',text)
 assert ast.literal_eval(re.search(r'NONCUSP counts deg 1\.\.9: (\[.*?\])',text)[1])==nc
 assert ast.literal_eval(re.search(r'REPS up to diamond: (\[.*?\])',text)[1])==reps
 assert 'classgroup=[ 8766, 0 ] rational_cusp_image=[ 4383 ] quotient=[ 2, 0 ]'in text
 assert 'complement_counts=[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]'in text
assert re.search(r'^(SIEVE_DONE|SIEVE_RESUME_DONE FromDegree=9 MaxPart=6)$',logs[-1],re.M) and 'total_survivors=0'in logs[-1]
ind=(H/'independent_counts_D9.log').read_text();assert 'DEGREE9_FAST_COUNTS_DONE complete=true'in ind
row=re.search(r'DEGREE9_FAST_COUNT m=1 n=38 p=3 genus=28 .*?closed_place_counts=(\[.*?\])',ind)
assert row and ast.literal_eval(row[1])==nc

def partitions(n,m=None):
 if n==0:yield ();return
 for a in range(min(n,m or n),0,-1):
  for tail in partitions(n-a,a):yield(a,)+tail
complete=[{tuple(ast.literal_eval(a)):(int(n),int(s))for a,n,s in re.findall(r'^RESULT type=(\[.*?\]) total=(\d+) survivors=(\d+)',t,re.M)}for t in logs]
empty=[{tuple(ast.literal_eval(a))for a in re.findall(r'^TYPE (\[.*?\]): no divisors',t,re.M)}for t in logs]
rows=[];zero=0
for d in range(1,10):
 for part in partitions(d):
  rest=Counter(part[1:]);expected=reps[part[0]-1]*prod(comb(nc[e-1]+k-1,k)if nc[e-1]else 0 for e,k in rest.items())
  if expected==0:assert any(part in a for a in empty);zero+=1;continue
  chosen=next((i for i,r in enumerate(complete)if r.get(part)==(expected,0)),None);assert chosen is not None,part
  cases=[int(i)for p,i,dead in re.findall(r'^SMART_CASE_RESULT type=(\[.*?\]) index=(\d+) killed=(\w+)',logs[chosen],re.M)if tuple(ast.literal_eval(p))==part and dead=='true']
  assert cases==list(range(1,expected+1)),(part,cases)
  rows.append({'partition':part,'count':expected,'completed_log':names[chosen]})
assert len(rows)==10 and zero==86 and sum(r['count']for r in rows)==438
control=(H/'N38_native_Hecke_control.log').read_text()
assert 'NONIC38_NATIVE_HECKE_EQUAL q=5 degree=54 principal=false'in control
assert 'NONIC38_NATIVE_CONTROL_COMPLETE'in control and not re.search('Runtime error|Assertion failed',control)
report={'degree':9,'excluded_group':[1,38],'coverage_complete':True,'noncuspidal_counts':nc,'diamond_representative_counts':reps,'partition_types':96,'empty_types':86,'nonempty_types':rows,'tested_representatives':438,'all_survivors':0,'untouched_native_Hecke_control':True,'audited_logs':names}
(H.parent/'data'/'N38_D9_audit.json').write_text(json.dumps(report,indent=2)+'\n');print('NONIC_C38_COVERAGE_PASS types=96 nonempty=10 tested=438 survivors=0')
