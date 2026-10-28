#!/usr/bin/env python3
"""Read-only audit of the complete nonic C35 exclusion."""
from pathlib import Path
import gzip,hashlib,json,re,subprocess,sys
P=Path(__file__).resolve().parent;ROOT=P.parents[1];D=ROOT/'C35/data/presentation_v2';LG=ROOT/'C35/logs/presentation_v2'
A=json.loads((D/'native_proof_audit.json').read_text());assert A['native_proven'] and A['factor_basis_bound']==A['unconditional_generation_bound_q3_g25']==9
assert 'bound=9' in (LG/'check_generation_bound.log').read_text()
B=json.loads((D/'portable_input_binding.json').read_text());assert hashlib.sha256(gzip.decompress((D/'proven_pic_presentation.json.gz').read_bytes())).hexdigest()==B['native_export_sha256']
assert hashlib.sha256((D/'verification_inputs.m').read_bytes()).hexdigest()==B['literal_input_sha256']
for f in sorted((ROOT/'C35/logs').rglob('*.log')):
 L=f.read_text();assert not re.search(r'Runtime error|User error|Assertion failed|Segmentation fault',L),f.name
C=json.loads((D/'cusp35.json').read_text());assert C['complete'] and C['global_and_F3_orbit_partitions_equal'] and C['all_F3_cusps_lift_to_same_global_degree']
assert C['global_signature']==C['F3_signature']=={'1':12,'4':3,'6':2,'12':1}
assert 'Gon_Q(X1(35)) = 12.' in (ROOT/'common/data/C35/gonality').read_text()
assert 'NONIC35_GLOBAL_OCTIC_VERIFIED degree8=true exact_P35=true F3_inert=true' in (ROOT/'C35/logs/verify_global_octic.log').read_text()
L=(LG/'restore_classify_dual_FB.log').read_text();assert len(re.findall('NONIC35_FB_RESTORE_PROGRESS',L))==14
assert 'NONIC35_FB_PRIMES_RESTORED places=3592 all_degrees_le9=true universal_generation_bound=9' in L
assert 'NONIC35_PRIMARY_RATIONAL_IMAGES full2_order=16384 full5_order=25 cusp2_order=8 cusp5_order=5 quotient_index=10240 exact=true' in L
assert 'NONIC35_PRIMARY_CLASSIFICATION_DONE tests=5961 unexpected=0 known_octic_places=6 complete=true' in L
assert 'native_proven=true' in (LG/'build.log').read_text()
for src,marker in [('replay_dual35.py','NONIC35_PRIMARY_INTEGER_REPLAY'),('replay_full_classes35.py','NONIC35_INDEPENDENT_CLASS_REPLAY')]:
 R=subprocess.run([sys.executable,str(P/src)],capture_output=True,text=True,check=True)
 assert marker in R.stdout and '"supports_exclusion": true' in R.stdout
print('NONIC_C35_READONLY_CERTIFICATE_PASS degree=9 native_proven=true generation_bound=9 FB_places=3592 native_relation_rows=4519 primary_image_orders=16384,25 rational_primary_orders=8,5 support_tests=5961 only_six_known_octic_places=true complete=true')
