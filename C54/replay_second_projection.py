#!/usr/bin/env python3
from pathlib import Path
import ast,json,sys
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent/'data'
sys.path.insert(0,str(P.parent));from replay_complete_classes import load,check_context,replay_rows
D=load('projection2_v2/results.py');assert D['complete'];check_context(D)
assert len(D['rows'])==13
labels=[]
for row in D['rows']:
 label,degrees,image_degrees,aa,member,passes,im,pre=row
 assert len(aa)==4 and aa[3]==0 and aa[2]%6518!=0
 assert member==passes==0 and im==pre==[];assert sum(degrees)==9
 assert all(d%e==0 for d,e in zip(degrees,image_degrees))
 labels.append(tuple(label))
assert labels==[(0,1,1,1),(0,2,1,1),(0,2,1,2),(0,2,1,3),(0,3,1,1),(0,3,2,1),(0,3,3,1),(1,6913,1),(1,7101,1),(1,6206,1),(1,5164,1),(1,5164,2),(1,5164,3)]
M=load('mixed_inert23/results.py');assert M['complete'];R=replay_rows(M);assert R['tests']==30 and R['members']==R['passes']==0
S=load('mixed_inert23/second_results.py');assert S=={'rows':[],'complete':True}
C=json.loads((P/'complete_classes_integer_replay.json').read_text());assert C['pure9_survivors']==[5164,6206,6913,7101]
def magma_values(name):
 out={}
 for line in (P/name).read_text().splitlines():
  if ':=' in line:
   key,value=line.split(':=',1);out[key]=ast.literal_eval(value.rstrip(';'))
 return out
I=magma_values('projection2_v2/inputs.m')
assert I['LowerDescriptors']==C['lower']['survivors']
assert set(I['Pure9Indices'])==set(C['pure9_survivors'])
L=load('lower_testclasses.py');first={tuple(row[0]):row for row in L['tests'] if row[3]}
for desc,pre in zip(I['LowerDescriptors'],I['LowerPreimages']):assert first[tuple(desc)][5]==pre
pure={row[0]:row for sh in range(4) for row in load(f'filtered9_classes_shard{sh}.py')['tests'] if row[3]}
for idx,pre in zip(I['Pure9Indices'],I['Pure9Preimages']):assert pure[idx][5]==pre
report={'mixed_inert23_tests':30,'mixed_inert23_members':0,'second_projection_tests':13,'all_second_classes_outside_H5':True,'all_first_survivors_exactly_bound_to_native_inputs':True,'complete':True}
if '--write-json' in sys.argv:(P/'second_projection_integer_replay.json').write_text(json.dumps(report,indent=2)+'\n')
print('NONIC54_SECOND_COMPLETE_INTEGER_REPLAY '+json.dumps(report))
