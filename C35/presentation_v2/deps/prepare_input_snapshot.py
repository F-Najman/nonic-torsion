#!/usr/bin/env python3
"""Translate a portable candidate workspace into literal Magma proof inputs.
No class-group claim is inferred from the candidate's native_proven flag.
"""
from pathlib import Path
import hashlib,json,sys
P=Path(sys.argv[1]).resolve();D=json.loads(P.read_text());assert D['checkpoint_complete']
assert D['status']=='candidate_presentation_only'
def m(a):
 if isinstance(a,bool):return str(a).lower()
 if isinstance(a,int):return str(a)
 return '['+','.join(map(m,a))+']'
fb=D['factor_basis'];out=P.with_name('verification_inputs.m')
s='/* Literal portable candidate input; not yet a proved Picard presentation. */\n'
for k,key in [('N','N'),('p','p'),('ExpectedGenus','genus'),('KnownOrder','known_Pic_order'),('FieldPolynomialData','field_polynomial'),('ModularCoordinatesData','modular_coordinates'),('FrobeniusData','frobenius_coefficients'),('CandidateInvariants','candidate_invariants'),('RawRelations','relations'),('DegreeOneVector','degree_one_vector'),('NaturalPoleVector','natural_pole_vector')]:s+=k+':='+m(D[key])+';\n'
s+='FBDegrees:='+m([a['degree'] for a in fb])+';\n'
s+='FBFinite:='+m([a['finite'] for a in fb])+';\n'
s+='FBBases:=[*'+','.join(m(a['ideal_basis']) for a in fb)+'*];\n'
out.write_text(s)
driver=out.with_name('verify_driver.m')
driver.write_text('OutDir:=\"'+str(P.parent)+'\";\nload \"'+str(out)+'\";\nload \"'+str(Path(__file__).resolve().parent/'verify.m')+'\";\n')
status={'candidate_file':P.name,'candidate_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'candidate_bytes':P.stat().st_size,'input_file':out.name,'input_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'claim':'candidate_only_until_explicit_relations_and_prime_lines_verified'}
P.with_name('verification_input_status.json').write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status))
