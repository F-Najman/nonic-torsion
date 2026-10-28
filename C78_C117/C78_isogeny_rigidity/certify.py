"""Read-only complete C78 exclusion certificate.

The new P9 Magma producer and its currently paused state are NOT premises.
The proof uses ordinary reduction and an odd-field isogeny cycle, after the
accepted degree-independent X139 allbad/lower-support/dyadic arguments.
All independent arithmetic replays are light and run by default.
"""
from pathlib import Path
import json,hashlib,subprocess,os,sys,time,re,math
P=Path(__file__).resolve().parent;D=P.parents[1];m=json.loads((D/'C78_C117/data/C78_isogeny_rigidity/manifest.json').read_text());assert m['degree']==9 and m['complete_exclusion']and m['excluded_groups']==[[1,78]]and not m['new_P9_Magma_filter_used']and not m['CM_classification_used']and not m['known_nonic39_model_used']
S=D/'C78_C117/C78_X139';O=D/'C78_C117/data/degree8';old=D/'C78_C117/degree8';om=json.loads((O/'quotient_manifest.json').read_text());assert hashlib.sha256((O/'quotient_manifest.json').read_bytes()).hexdigest()=='809e63dc15648d1e3af406d6503e5afee7f6b3aea9407987f8cda702d525bb71'
assert om['complete']and om['source_rational_Jacobian_invariants']==[7,31122,3236688]and om['rational_gonality']==14 and om['target_rational_cusp_2primary_image']==32 and om['good_nonprincipality_witnesses']==5 and om['allbad_cusps_lift_Q']
s=(D/'C78_C117/logs/degree8/hecke_good_39.log').read_text();assert'C39_D8_GOOD_DONE'in s and'unknown=0 complete=true'in s and not any(x in s for x in['Runtime error','User error','Assertion failed']);w=re.findall(r'C39_D8_GOOD_WITNESS degree=(\d+) index=(\d+)/(\d+) q=(\d+) principal=(true|false)',s);assert len(w)==5 and all(q=='5'and b=='false'for _,_,_,q,b in w)
assert'HeckeOperatorCommonZeros(X,q,DD)-q*DiamondOperator(X,q,DD)-DD'in(old/'hecke_good_39.m').read_text()
R=D/'C78_C117/data/C78_X139'
c=json.loads((R/'census39_F512.json').read_text());assert c['complete']and sorted(v['order']for v in c['all_supersingular_four_models_directly_counted'])==[481,513,513,545]and all(v['order']%39 for v in c['all_supersingular_four_models_directly_counted'])
assert 14>9 and all(r%2==1 and math.isqrt(2**r)**2!=2**r for r in[1,3,9])
commands=[([sys.executable,str(S/'replay_dyadic_injection39.py')],'C39_C4_INDEPENDENT_REPLAY_DONE'),([sys.executable,str(S/'census39_F512.py'),'--readonly'],'C39_F512_INDEPENDENT_COMPLETE_CENSUS_PASS')]
env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
for cmd,marker in commands:
 t=time.monotonic();r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,env=env);assert r.returncode==0 and marker in r.stdout,(cmd,r.returncode,r.stdout);print('C78_LIGHT_READ_ONLY_REPLAY_PASS',Path(cmd[1]).name,round(time.monotonic()-t,3),flush=True)
print('C78_D9_COMPLETE_READ_ONLY_CERTIFICATE_PASS excluded=C78 new_P9_Magma_filter_used=false CM_classification_used=false known_nonic39_model_used=false files='+str(len(m['files'])),flush=True)
