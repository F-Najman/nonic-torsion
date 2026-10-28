from pathlib import Path
import subprocess,resource,time,json,sys
P=Path(__file__).resolve().parent;ROOT=P.parent;L=ROOT/'C2x42/logs'
def limits():resource.setrlimit(resource.RLIMIT_CPU,(14400,14430))
records=[]
for d in range(4,10):
 start=time.monotonic();log=L/('enumerate_d'+str(d)+'.log')
 with log.open('w') as f:r=subprocess.run(['magma','-n','-b','D:='+str(d),'Limit:=0','Tag:=full'+str(d),str(P/'enumerate_direct.m')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits)
 s=log.read_text();assert r.returncode==0 and 'M2N42_D9_DIRECT_DONE' in s and not any(x in s for x in ['Runtime error','User error','syntax error','Assertion failed']),s[-2000:]
 records.append(dict(degree=d,actual_exit_code=r.returncode,elapsed_seconds=time.monotonic()-start));(P/'enumeration_completion.json').write_text(json.dumps(dict(completed_degrees=[r['degree'] for r in records],complete=d==9,records=records),indent=2)+'\n')
 print('C242_D9_DIRECT_NATIVE_PASS degree='+str(d),flush=True)
