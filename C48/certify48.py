#!/usr/bin/env python3
"""Read-only computational certificate for C48 over degree9 fields."""
from pathlib import Path
import os,sys,json,hashlib,subprocess,runpy
sys.dont_write_bytecode=True
from certificate_common import literal,magma
P=Path(__file__).resolve().parent
# Bind the new portable presentation to actual archived, proved finite data.
G=literal(P/'data/degree8/clean_run/factorbase_geometry.py');B=literal(P/'data/degree8/clean_run_v2/function_field_bindings.py');A=literal(P/'data/degree8/clean_run/target_2_24_arithmetic.py');V=magma(P/'data/portable_inputs.m')
expected={'PortableFieldEquation':G['field_equation'],'PortableFlags':G['finite_flags'],'PortableDegrees':G['factorbase_degrees'],'PortableIdealBases':G['factorbase_ideal_bases'],'PortableClasses':G['factorbase_classes'],'PortableRelations':G['relations'],'PortableCoordinates':B['target_modular_coordinates'],'PortableBaseSignature':B['base_cusp_signature'],'PortableA':[[801,1]],'PortableRationalCusps':A['rational_cusp_generators'],'PortableCuspClasses':A['finite_cusp_classes'],'PortableCuspDegrees':A['finite_cusp_degrees']}
assert V==expected
# Bind character inputs and all actual class-query RR witnesses.
C=literal(P/'data/finite_characters_v2/function_data.py');F=magma(P/'data/character_filter9/functions.m');assert F=={'FilterFunctions':C['functions'],'FilterBaseValues':C['base_values'],'FilterCharacterRows':C['character_rows'],'FilterCombinations':C['annihilating_combinations']}
R=magma(P/'data/parameter9/records.m');assert R['Complete'] and len(R['Records'])==13513 and sum(v[4] for v in R['Records'])==217520 and sum(v[3] for v in R['Records'])==216208
saved=[]
for sh in [0,1]:
 a=literal(P/f'data/class9_shard{sh}/results.py');w=literal(P/f'logs/class9_shard{sh}/identities.log.witnesses.py')['witnesses'];assert len(a['tests'])==len(w)==348 and a['complete']
 for r,v in zip(a['tests'],w):assert v[1]==0 and v[2]==r[1] and v[5]==[[801,1]];saved.append([r[0],r[1],v])
saved.sort(key=lambda r:r[0]);assert magma(P/'data/class9_witness_verifier_v2/inputs.m')['ReplayCases']==saved
for b in json.loads((P/'data/class9_witness_verifier_v2/input_bindings.json').read_text())['sources']:assert hashlib.sha256((P/('logs' if 'witnesses' in b['path'] else 'data')/b['path']).read_bytes()).hexdigest()==b['sha256']
env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
for prog in ['replay_D9_integer.py','replay_source_Hecke_coverage.py','replay_character9_sage.py']:
 r=subprocess.run([sys.executable,str(P/prog)],capture_output=True,text=True,env=env);assert r.returncode==0,(prog,r.stdout,r.stderr);print(r.stdout.strip(),flush=True)
print('NONIC48_CERTIFICATE_PASS complete_exclusion=true remaining_source_cases=0',flush=True)
