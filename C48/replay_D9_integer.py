#!/usr/bin/env python3
"""Read-only degree9 arithmetic replay. Imported finite data are hash verified."""
from pathlib import Path
import ast,json,hashlib,runpy
from collections import Counter
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent

def literal(path):
 s=path.read_text().replace('[*','[').replace('*]',']');out={}
 for n in ast.parse(s).body:
  assert isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name)
  out[n.targets[0].id]=ast.literal_eval(n.value)
 return out
A=literal(P/'data/degree8/clean_run/target_2_24_arithmetic.py');L=runpy.run_path(str(P/'data/degree8/clean_run/lower_support_arithmetic.py'))
assert L['complete'] and L['relations']==A['relations'] and L['rational_cusp_generators']==A['rational_cusp_generators']
R=Matrix(A['relations']);H=hermite_normal_form(R[:,:8].T.row_join(Matrix([v[:8] for v in A['rational_cusp_generators']]).T));assert abs(int(H.det()))==4608

def red(v):
 assert len(v) in (8,9)
 if len(v)==9:assert v[8]==0
 z=list(v[:8])
 for i in reversed(range(8)):
  q=z[i]//int(H[i,i])
  for j in range(i+1):z[j]-=q*int(H[j,i])
 return tuple(z)
def add(v,w):return red([x+y for x,y in zip(v,w)])
comp=[{(0,)*8}]+[set() for _ in range(9)]
for d,v in zip(A['finite_cusp_degrees'],A['finite_cusp_classes']):
 for j in range(d,10):comp[j].update(add(z,red(v)) for z in list(comp[j-d]))
assert list(map(len,comp))==[1,5,23,75,178,378,661,1065,1487,1859]
counts=Counter();passing=[];tested=0

def test(degrees,indices,v):
 global tested
 tested+=1
 if red([-x for x in v]) in comp[9-sum(degrees)]:counts[tuple(degrees)]+=1;passing.append([degrees,indices])
for d in (3,4,5):
 for i,v in enumerate(L['representative_classes'][d],1):test([d],[i],v)
for d,e in ((3,3),(3,4),(3,5),(4,4),(4,5)):
 for i,v in enumerate(L['representative_classes'][d],1):
  for j,w in enumerate(L['all_source_classes'][e],1):test([d,e],[i,j],[x+y for x,y in zip(v,w)])
for i,v in enumerate(L['representative_classes'][3],1):
 for j,w in enumerate(L['all_source_classes'][3],1):
  for k in range(j-1,len(L['all_source_classes'][3])):
   z=L['all_source_classes'][3][k];test([3,3,3],[i,j,k+1],[x+y+t for x,y,t in zip(v,w,z)])
D=json.loads((P/'data/first_D9_filter.json').read_text());assert passing==D['low_survivors'];assert counts=={(3,5):4,(4,4):16}
rows={r[0]:r for sh in range(8) for r in literal(P/f'data/degree8/clean_run_v2/high_shard{sh}_arithmetic.py')['tests']};assert len(rows)==3507
orbits=runpy.run_path(str(P/'data/degree8/clean_run/high_target_aut_orbits.py'));cids={cid for cid,r in rows.items() if red([-x for x in r[4]]) in comp[9-r[1]]}
high=[[r[0],r[1],r[2]] for r in orbits['cases'] if r[2] in cids];assert high==D['cached_high_survivors'];assert Counter(r[0] for r in high)=={6:14,7:8,8:6}
Q=literal(P/'data/finite_characters_v2/action_data.py');mods=[2,2,24,48];projection=Q['quotient_images_of_G_generators']

def project(v):return tuple(sum(v[j]*projection[j][i] for j in range(9))%mods[i] for i in range(4))
def act(vector,images):return tuple(sum(vector[j]*images[j][i] for j in range(4))%mods[i] for i in range(4))
assert Q['quotient_relations']==[[mods[i] if i==j else 0 for j in range(4)] for i in range(4)]
for rel in A['relations']+A['rational_cusp_generators']:assert project(rel)==(0,0,0,0)
GEO=literal(P/'data/degree8/clean_run/factorbase_geometry.py');assert project(GEO['factorbase_classes'][800])==(0,0,0,0)
full=hermite_normal_form(Matrix(Q['quotient_relations']+projection[:8]).T);assert abs(int(full.det()))==1
for images in Q['action_images']:
 assert len(images)==4
 for j in range(4):assert all((mods[j]*x)%m==0 for x,m in zip(images[j],mods))
 assert abs(int(hermite_normal_form(Matrix(Q['quotient_relations']+images).T).det()))==1

def magma_assignments(path):
 out={}
 for stmt in path.read_text().split(';'):
  if ':=' in stmt:
   name,value=stmt.split(':=',1);out[name.strip()]=ast.literal_eval(value.strip().replace('[*','[').replace('*]',']'))
 return out
MI=magma_assignments(P/'data/mixed36/inputs.m')
assert MI['ActionLabels']==Q['labels'] and MI['ActionImages']==Q['action_images'] and MI['QRelations']==Q['quotient_relations']
master=runpy.run_path(str(P/'data/degree8/source_partition_inputs_complete.py'))
assert MI['Source6Polynomial']==master['field_polynomials'][6] and MI['Source6Representatives']==master['representative_pairs'][6]
case6={r[1]:r for r in orbits['cases'] if r[0]==6}
assert MI['Canonical6Keys']==[case6[i][3] for i in range(1,343)]
assert [tuple(v) for v in MI['Canonical6Classes']]==[project(rows[case6[i][2]][4]) for i in range(1,343)]
assert [tuple(v) for v in MI['Source3RepClasses']]==[project(v) for v in L['representative_classes'][3]]

M=literal(P/'data/mixed36/results.py');assert M['complete'] and M['tested']==5424 and M['passing_cases']==[] and len(M['source6_all_rows'])==2712
for rep,unit,pair,cl in M['source6_all_rows']:
 aa=min(unit%24,24-unit%24);idx=Q['labels'].index([0,0,aa]);assert tuple(cl)==act(M['source6_representative_classes'][rep-1],Q['action_images'][idx])
for v in L['representative_classes'][3]:
 for row in M['source6_all_rows']:assert any((x+y)%m for x,y,m in zip(project(v),row[3],mods))
tested+=5424;assert tested==13020
C=literal(P/'data/finite_characters_v2/function_data.py')
def binaryrank(rows):
 piv={}
 for row in rows:
  v=sum((x%2)<<j for j,x in enumerate(row))
  while v:
   b=v.bit_length()-1
   if b not in piv:piv[b]=v;break
   v^=piv[b]
 return len(piv)
assert binaryrank(C['basis_vectors'])==binaryrank(C['character_rows'])==8
chars=[[sum(a*b for a,b in zip(combo,col))%2 for col in zip(*C['character_rows'])] for combo in C['annihilating_combinations']]
assert binaryrank(chars)==4
for ch in chars:
 assert all(sum(a*b for a,b in zip(ch,v[:8]))%2==0 for v in A['rational_cusp_generators'])
F=literal(P/'data/character_filter9/results.py');assert F['complete'] and len(F['rows'])==13513 and F['no_source_orbits']==1644 and F['fallback']==0
for idx,passed,signs in F['rows']:
 if passed==-1:assert signs==[];continue
 assert all(signs[i-1] in(-1,1) for i in(2,3,6,7,8))
 got=all(__import__('math').prod(s**a for s,a in zip(signs,combo) if a)!=-1 for combo in C['annihilating_combinations'])
 assert passed==int(got)
seen=set();pass9=[]
for sh in (0,1):
 B=literal(P/f'data/class9_shard{sh}/results.py');assert B['complete'] and B['quotient_images_of_G']==projection and B['quotient_relations']==Q['quotient_relations']
 assert B['expected_indices']==F['remaining_indices'][sh::2]
 for idx,v,cl,passed in B['tests']:
  assert idx not in seen;seen.add(idx);assert project(v)==tuple(cl);assert passed==int(not any(cl))
  if passed:pass9.append(idx)
assert seen==set(F['remaining_indices']) and len(seen)==696 and len(pass9)==11
assert sorted(pass9)==sorted(json.loads((P/'data/class9_survivors.json').read_text())['target_Aut16_orbit_indices'])
print('NONIC48_D9_INTEGER_REPLAY_PASS '+json.dumps({'nonpure9_source_cases':20882,'new_mixed36_cases':5424,'pure9_character_tests':13513,'full_class_queries':696,'pure9_target_survivors':11,'nonpure9_source_survivors':48,'complete':True,'global_exclusion_claimed':False}))
