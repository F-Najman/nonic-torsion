#!/usr/bin/env python3
"""Read-only complete order93 exclusion replay; never launches Magma."""
from pathlib import Path
import os,sys,json,hashlib,subprocess,resource,ast,re,time,math
K=Path(__file__).resolve().parent;H=K.parent/'data/N93_source_units';D=K.parents[1];ROOT=D;t=time.monotonic()
m=json.loads((H/'manifest93.json').read_text());assert m['complete_exclusion']and m['group']==[1,93]
s=(D/'C85_C91_C93/data/character_probes/N93_winding_data.m').read_text()
def lit(n):return ast.literal_eval(re.search(re.escape(n)+r':=(.*?);',s,re.S)[1].replace('[*','[').replace('*]',']').replace('<','(').replace('>',')'))
f=lit('ConstituentWindingFlags');z=[v for v in f if v[2]==0];assert len(z)==1 and z[0][0]==2 and z[0][1]==93 and z[0][3]==1
assert lit('FiniteProjectiveDiamondOrder')==30 and lit('FiniteDiamondUnits')==[a for a in range(1,94)if math.gcd(a,93)==1]
assert pow(4,5,93)==1 and all(pow(4,j,93)not in[1,92]for j in range(1,5))
q=json.loads((H/'canonical_R1_sieve_verified.json').read_text());assert q['status']=='PASS'and q['cases']==120 and not q['unresolved']and all(v['rank']==261 and v['positive_degree']==108 and v['columns']==627 for v in q['results'])
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
def limits():
 os.nice(10);resource.setrlimit(resource.RLIMIT_CPU,(600,610));resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
checks=[['verify_poles93_readonly.py'],['verify_stage1_readonly_v3.py'],['verify_stage2_readonly_v2.py'],['build_R1_product_tensor.py','--readonly'],['canonical_R1_sieve93.py','--readonly'],['verify_R1_ranks_m4ri.py']]
for args in checks:
 r=subprocess.run([sys.executable,'-B',str(K/args[0])]+args[1:],cwd=ROOT,env=env,preexec_fn=limits,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=900)
 if r.returncode:print(r.stdout);raise RuntimeError((args[0],r.returncode))
 print('READONLY_CHECK_PASS',args[0],flush=True)
print('NONIC93_COMPLETE_READONLY_CERTIFICATE_PASS',json.dumps({'files':len(m['files']),'native_and_independent_successful_jobs':len(m['successful_jobs']),'allbad_gonality':'>36','source_points':10,'actual_isogenies_independently_rebuilt':6,'binary_matrices':120,'rank':261,'elapsed_seconds':time.monotonic()-t}))
