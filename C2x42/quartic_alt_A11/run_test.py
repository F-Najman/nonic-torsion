from pathlib import Path
import subprocess,resource,os,time,json,sys,signal
P=Path(__file__).resolve().parent;ROOT=P.parents[1];L=ROOT/'C2x42/logs/quartic_alt_A11'
def limits():
 resource.setrlimit(resource.RLIMIT_CPU,(1770,1800));resource.setrlimit(resource.RLIMIT_AS,(10*10**9,10*10**9));os.nice(10)
t=time.monotonic();timed=False
with (L/'test.log').open('w') as log:
 proc=subprocess.Popen(['magma','-n','-b',str(P/'test.m')],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits,start_new_session=True)
 (P/'launch.json').write_text(json.dumps(dict(supervisor_pid=os.getpid(),launcher_pid=proc.pid,cpu_soft=1770,cpu_hard=1800,memory_budget=8*10**9,address_space_limit=10*10**9,wall_limit=7200,nice=10),indent=2)+'\n')
 time.sleep(.1)
 q=Path('/proc')/str(proc.pid)
 if proc.poll() is None and q.exists():
  st=(q/'stat').read_text().split();j=dict(native_pid=proc.pid,start_ticks=int(st[21]),command=(q/'cmdline').read_bytes().replace(b'\x00',b' ').decode(),memory_budget=8*10**9,cpu_hard=1800)
  (P/'native_identity.json').write_text(json.dumps(j,indent=2)+'\n')
 try:rc=proc.wait(timeout=7200)
 except subprocess.TimeoutExpired:
  timed=True;os.killpg(proc.pid,signal.SIGTERM)
  try:rc=proc.wait(timeout=15)
  except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);rc=proc.wait()
s=(L/'test.log').read_text();out=dict(actual_exit_code=rc,wall_seconds=time.monotonic()-t,wall_limit=7200,wall_timeout=timed,cpu_hard=1800,memory_budget=8*10**9,complete=rc==0 and not timed and 'C242_QUARTIC_ALT_A11_DONE' in s)
(P/'completion.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True);sys.exit(0 if out['complete'] else 1)
