from pathlib import Path
import subprocess,resource,time,json,sys
P=Path(__file__).resolve().parent;ROOT=P.parents[0];L=ROOT/'C2x42/logs'
def limits():resource.setrlimit(resource.RLIMIT_CPU,(14400,14430))
t=time.monotonic()
with (L/'filter_quotient14_linear.log').open('w') as f:r=subprocess.run(['magma','-n','-b',str(P/'filter_quotient14_linear.m')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits)
s=(L/'filter_quotient14_linear.log').read_text();a=dict(actual_exit_code=r.returncode,elapsed_seconds=time.monotonic()-t,complete=r.returncode==0 and 'C242_D9_QUOTIENT14_LINEAR_DONE' in s)
(P/'filter_quotient14_linear_completion.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a));sys.exit(0 if a['complete'] else 1)
