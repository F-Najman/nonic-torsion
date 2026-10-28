#!/usr/bin/env python3
"""Read-only certificate for the degree9 C54 exclusion."""
from pathlib import Path
import os,json,hashlib,subprocess,sys
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;LG=P/'logs';DA=P/'data'
for name,marker in [('projection2_v2','NONIC54_SECOND_SURVIVOR_TEST_DONE'),('mixed_inert23','NONIC54_MIXED_INERT23_DONE'),('positive7_control','NONIC54_RETAINED7_POSITIVE_DONE')]:
 log=(LG/name/'run.log').read_text();assert marker in log
for f in sorted(LG.rglob('*.log')):
 s=f.read_text();assert all(z not in s for z in ['Runtime error','User error','Assertion failed','Segmentation fault']),f.name
# The inherited degree8 class packet is used only for finite H5 nonmembership;
# every positive representative is freshly verified in this degree9 directory.
old=json.loads((DA/'degree8/old_root_manifest.json').read_text())
assert old['degree']==8 and old['group']==[1,54] and old['complete_exclusion']
for sh in range(4):
 for ext in ['py','log']:
  original=f'large_levels/C54_good8_classes_shard{sh}.{ext}' if ext=='py' else f'large_levels/C54_good8_class_shard{sh}.{ext}'
  assert hashlib.sha256(((DA if ext=='py' else LG)/f'degree8/old_good8_classes_shard{sh}.{ext}').read_bytes()).hexdigest()==old['files'][original]
env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
for prog,args in [('replay_class_arithmetic.py',['preflight']),('replay_quadratic_filter.py',[]),('replay_complete_classes.py',[]),('replay_second_projection.py',[])]:
 r=subprocess.run([sys.executable,str(P/prog),*args],env=env,capture_output=True,text=True)
 assert r.returncode==0,(prog,r.stdout,r.stderr)
 print(r.stdout.strip())
F=json.loads((DA/'final_checkpoint54.json').read_text());assert F['final_survivors']==0 and F['second_projection_H5_members']==0
print('NONIC54_CERTIFICATE_PASS complete_exclusion=true final_survivors=0')
