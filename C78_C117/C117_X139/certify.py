"""Complete read-only C117 exclusion certificate; no new Magma process.

The still-paused general X139 nonic Hecke producer is NOT a premise.
Only the accepted degree-independent lower-degree/dyadic certificate,
the independent complete F512 census, and the known39 positive models
are used. Every independent arithmetic replay is light and runs by default.
"""
from pathlib import Path
import json,hashlib,subprocess,os,sys,time,re
P=Path(__file__).resolve().parent;D=P.parents[1];m=json.loads((D/'C78_C117/data/C117_X139/manifest.json').read_text());assert m['degree']==9 and m['complete_exclusion']and m['excluded_groups']==[[1,117]]and not m['new_P9_Magma_filter_used']
S=D/'C78_C117/C78_X139';O=D/'C78_C117/data/degree8';old=D/'C78_C117/degree8';om=json.loads((O/'quotient_manifest.json').read_text());assert hashlib.sha256((O/'quotient_manifest.json').read_bytes()).hexdigest()=='809e63dc15648d1e3af406d6503e5afee7f6b3aea9407987f8cda702d525bb71'
assert om['complete']and om['source_rational_Jacobian_invariants']==[7,31122,3236688]and om['rational_gonality']==14 and om['target_rational_cusp_2primary_image']==32 and om['good_nonprincipality_witnesses']==5 and om['allbad_cusps_lift_Q']
s=(D/'C78_C117/logs/degree8/hecke_good_39.log').read_text();assert'C39_D8_GOOD_DONE'in s and'unknown=0 complete=true'in s and not any(x in s for x in['Runtime error','User error','Assertion failed']);w=re.findall(r'C39_D8_GOOD_WITNESS degree=(\d+) index=(\d+)/(\d+) q=(\d+) principal=(true|false)',s);assert len(w)==5 and all(q=='5'and b=='false'for _,_,_,q,b in w)
assert'HeckeOperatorCommonZeros(X,q,DD)-q*DiamondOperator(X,q,DD)-DD'in(old/'hecke_good_39.m').read_text()
R=D/'C78_C117/data/C78_X139'
c=json.loads((R/'census39_F512.json').read_text());assert c['complete']and c['all511_ordinary_j_invariants_and_both_twists_checked']and c['all39_selected_orders_direct_character_sum_checked']and c['j_degree3_order468_curves']==3 and c['j_degree9_order546_curves']==36 and sorted(v['order']for v in c['all_supersingular_four_models_directly_counted'])==[481,513,513,545]
assert[n for n in range(468,559)if n%117==0]==[468]
K=D/'C78_C117/data/C78_known39';assert hashlib.sha256((K/'manifest39.json').read_bytes()).hexdigest()=='921fd3c267f053962ea1fb90c8936ac39be782114823b116c03c5c69e4e9ab05';j=json.loads((K/'positive_divisors39.json').read_text());assert j['global_closed_nonic_divisors']==4 and len(j['known_positive_divisors'])==4 and len({tuple(v['frobenius_canonical_Tate_bc_key'])for v in j['known_positive_divisors']})==4 and all(v['global_trace_degree']==9 and v['reduction_code'][0]==9 for v in j['known_positive_divisors'])
jj=json.loads((K/'independent39.json').read_text());assert jj['complete']and jj['full_global_torsion_order']==39 and jj['F29_order']==39 and jj['F512_trace_count']==468 and jj['good_unramified_prime']==29
assert 4*9==3*12
commands=[([sys.executable,str(S/'replay_dyadic_injection39.py')],'C39_C4_INDEPENDENT_REPLAY_DONE'),([sys.executable,str(S/'census39_F512.py'),'--readonly'],'C39_F512_INDEPENDENT_COMPLETE_CENSUS_PASS'),([sys.executable,str(D/'C78_C117/C78_known39/certify39.py')],'NONIC78_KNOWN39_READONLY_CERTIFICATE_PASS'),([sys.executable,str(D/'C78_C117/C117_known39_binding/certify117.py')]+m['binding_checker_args'],m['binding_checker_marker'])]
env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
for cmd,marker in commands:
 t=time.monotonic();r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,env=env);assert r.returncode==0 and marker in r.stdout,(cmd,r.returncode,r.stdout);print('C117_LIGHT_READ_ONLY_REPLAY_PASS',Path(cmd[1]).name,round(time.monotonic()-t,3),flush=True)
print('C117_D9_COMPLETE_READ_ONLY_CERTIFICATE_PASS excluded=C117 no_new_P9_Magma_filter_used=true four_known9_lifts=true files='+str(len(m['files'])),flush=True)
