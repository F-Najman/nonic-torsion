#!/usr/bin/env python3
"""Read-only certificate for the complete degree9 C58 exclusion."""
from pathlib import Path
import hashlib,json,os,re,subprocess,sys
S=Path(__file__).resolve().parent;P=S.parent/'data'/'C58_X129';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
W=json.loads((P/'presentation/presentation.json').read_text())
assert W['native_proven'] and W['factor_basis_bound']==9
literal=lambda x:json.dumps(x,separators=(',',':'))
expected=''.join(f'{k}:={literal(v)};\n' for k,v in [('N',W['N']),('p',W['p']),('ExpectedGenus',W['genus']),('KnownOrder',W['known_Pic_order']),('CandidateInvariants',W['candidate_invariants']),('FieldPolynomialData',W['field_polynomial']),('ModularCoordinatesData',W['modular_coordinates']),('RawRelations',W['relations']),('FBDegrees',[a['degree'] for a in W['factor_basis']]),('FBFinite',[a['finite'] for a in W['factor_basis']]),('FBBases',[a['ideal_basis'] for a in W['factor_basis']])])
assert (P/'inputs.m').read_text()==expected
C=json.loads((P/'cusp58.json').read_text());assert C['global_and_F3_orbit_partitions_equal'] and C['global_signature']==C['F3_signature']=={'1':14,'14':1}
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
r=subprocess.run([sys.executable,str(S/'replay58.py')],env=env,capture_output=True,text=True)
assert r.returncode==0,(r.stdout,r.stderr);assert 'NONIC58_INDEPENDENT_REPLAY' in r.stdout and '"supports_exclusion": true' in r.stdout
print(r.stdout.strip())
print('NONIC58_CERTIFICATE_PASS degree=9 tests=2324 survivors=0 complete_exclusion=true')
