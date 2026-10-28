#!/usr/bin/env python3
"""Condition (pushclass) in degree 9 for the supports of types other than (9) and (3,6), from the class data of the
degree-8 filters; writes C48/data/first_D9_filter.json and the portable target presentation C48/data/portable_inputs.m.
Extracted for this repository from the inline Python command of the Codex session (transcript 2026-10-26T19:31:55Z)
that wrote both files; only the paths are adapted (the clean_run_v2 copies of target_2_24_arithmetic.py and
factorbase_geometry.py are byte-identical to the clean_run copies, so the first printed comparison is now trivial).
C48/replay_D9_integer.py recomputes first_D9_filter.json independently, and C48/certify48.py binds portable_inputs.m
to the degree-8 data."""
from pathlib import Path
import runpy,ast,collections,json
p=Path(__file__).resolve().parent;o=p/'data'/'degree8';lg=p/'logs'/'degree8'
def lit(f):return {n.targets[0].id:ast.literal_eval(n.value) for n in ast.parse(f.read_text()).body}
A=lit(o/'clean_run/target_2_24_arithmetic.py');B=lit(o/'clean_run/target_2_24_arithmetic.py');L=runpy.run_path(str(o/'clean_run/lower_support_arithmetic.py'))
print('presentations_equal',A['relations']==B['relations'],'H_equal',A['rational_cusp_generators']==B['rational_cusp_generators'],'low_equals_high',L['relations']==B['relations'] and L['rational_cusp_generators']==B['rational_cusp_generators'])
geom=lit(o/'clean_run/factorbase_geometry.py');binding=lit(o/'clean_run_v2/function_field_bindings.py');W=lit(lg/'clean_run_v2/high_shard0_identities.log.witnesses.py')['witnesses'];print('baseA',W[0][5])
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
R=Matrix(B['relations']);H=hermite_normal_form(R[:,:8].T.row_join(Matrix([v[:8] for v in B['rational_cusp_generators']]).T))
def red(v):
 z=list(v[:8])
 for i in reversed(range(8)):
  q=z[i]//int(H[i,i])
  for j in range(i+1):z[j]-=q*int(H[j,i])
 return tuple(z)
def add(v,w):return red([a+b for a,b in zip(v,w)])
comps=[{(0,)*8}]+[set() for _ in range(9)]
for deg,v in zip(B['finite_cusp_degrees'],B['finite_cusp_classes']):
 for j in range(deg,10):comps[j].update(add(z,red(v)) for z in list(comps[j-deg]))
print('D9_comps',list(map(len,comps)))
counts=collections.Counter();surv=[]
for d in [3,4,5]:
 for i,v in enumerate(L['representative_classes'][d],1):
  if red([-a for a in v]) in comps[9-d]:counts[str([d])]+=1;surv.append([[d],[i]])
for d,e in [(3,3),(3,4),(3,5),(4,4),(4,5)]:
 for i,v in enumerate(L['representative_classes'][d],1):
  for j,w in enumerate(L['all_source_classes'][e],1):
   if red([-a-b for a,b in zip(v,w)]) in comps[9-d-e]:counts[str([d,e])]+=1;surv.append([[d,e],[i,j]])
for i,v in enumerate(L['representative_classes'][3],1):
 for j,w in enumerate(L['all_source_classes'][3],1):
  for k in range(j-1,len(L['all_source_classes'][3])):
   z=L['all_source_classes'][3][k]
   if red([-a-b-c for a,b,c in zip(v,w,z)]) in comps[0]:counts[str([3,3,3])]+=1;surv.append([[3,3,3],[i,j,k+1]])
print('newlow_passcounts',dict(counts),'total',len(surv))
rows={r[0]:r for sh in range(8) for r in lit(o/f'clean_run_v2/high_shard{sh}_arithmetic.py')['tests']}
orb=runpy.run_path(str(o/'clean_run/high_target_aut_orbits.py'))
passcids={cid for cid,row in rows.items() if red([-a for a in row[4]]) in comps[9-row[1]]}
highsurv=[[r[0],r[1],r[2]] for r in orb['cases'] if r[2] in passcids]
print('newhigh_passcounts',collections.Counter(r[0] for r in highsurv),'total',len(highsurv))
(p/'data'/'first_D9_filter.json').write_text(json.dumps({'complement_counts':list(map(len,comps)),'low_survivors':surv,'low_counts':dict(counts),'cached_high_survivors':highsurv,'unknown_supports':[[3,6],[9]],'degree8_exclusion_transferred':False},indent=2)+'\n')
def magma(v):
 if isinstance(v,int):return str(v)
 if isinstance(v,list):return '['+','.join(map(str,v))+']' if all(isinstance(z,int) for z in v) else '[*'+','.join(magma(z) for z in v)+'*]'
 raise ValueError(type(v))
values={'PortableFieldEquation':geom['field_equation'],'PortableFlags':geom['finite_flags'],'PortableDegrees':geom['factorbase_degrees'],'PortableIdealBases':geom['factorbase_ideal_bases'],'PortableClasses':geom['factorbase_classes'],'PortableRelations':geom['relations'],'PortableCoordinates':binding['target_modular_coordinates'],'PortableA':W[0][5],'PortableBaseSignature':binding['base_cusp_signature'],'PortableRationalCusps':B['rational_cusp_generators'],'PortableCuspClasses':B['finite_cusp_classes'],'PortableCuspDegrees':B['finite_cusp_degrees']}
(p/'data'/'portable_inputs.m').write_text(''.join(k+':='+magma(v)+';\n' for k,v in values.items()))
print('inputs_written')
