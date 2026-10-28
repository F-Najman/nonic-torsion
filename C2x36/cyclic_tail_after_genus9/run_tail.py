from pathlib import Path
import subprocess,resource,time,json,sys
P=Path(__file__).resolve().parent;ROOT=P.parents[1];L=ROOT/'C2x36/logs/cyclic_tail_after_genus9'
def limits():resource.setrlimit(resource.RLIMIT_CPU,(14400,14430))
t=time.monotonic()
with (L/'tail.log').open('w') as f:r=subprocess.run(['magma','-n','-b',str(P/'tail.m')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits)
s=(L/'tail.log').read_text();a=dict(actual_exit_code=r.returncode,elapsed_seconds=time.monotonic()-t,complete=r.returncode==0 and 'C236_D9_REDUCED_CYCLIC36_TAIL_DONE' in s)
(P/'tail_completion.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a));sys.exit(0 if a['complete'] else 1)
