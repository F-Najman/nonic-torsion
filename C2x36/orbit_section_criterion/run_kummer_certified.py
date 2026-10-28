from pathlib import Path
import subprocess,resource,time,json,sys
P=Path(__file__).resolve().parent;ROOT=P.parents[1];L=ROOT/'C2x36/logs/orbit_section_criterion'
def limits():resource.setrlimit(resource.RLIMIT_CPU,(5400,5430))
t=time.monotonic()
with (L/'kummer_certified.log').open('w') as f:r=subprocess.run(['magma','-n','-b',str(P/'check_kummer_certified.m')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits)
s=(L/'kummer_certified.log').read_text();a=dict(actual_exit_code=r.returncode,elapsed_seconds=time.monotonic()-t,cpu_limit=5400,complete=r.returncode==0 and 'C236_D9_KUMMER_CERTIFIED_TRACE_ORBIT_SECTION_OBSTRUCTION_PASS' in s)
(P/'kummer_certified_completion.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a));sys.exit(0 if a['complete'] else 1)
