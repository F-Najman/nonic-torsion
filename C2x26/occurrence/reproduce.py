from pathlib import Path
import subprocess,sys,json,time,resource,argparse
P=Path(__file__).resolve().parent;ROOT=P.parents[1]
ap=argparse.ArgumentParser();ap.add_argument('output',type=Path);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
assert not any(c in str(out) for c in ['"','\n','\r'])
wrapper=out/'driver.m';wrapper.write_text('C226Out:='+json.dumps(str(out))+';\nload '+json.dumps(str(P/'finite_geometric_reproduction.m'))+';\n')
def limits():resource.setrlimit(resource.RLIMIT_CPU,(3600,3630))
start=time.monotonic()
with (out/'native.log').open('w') as f:r=subprocess.run(['magma','-b',str(wrapper)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits)
s=(out/'native.log').read_text();assert r.returncode==0 and 'C226_GEOMETRIC_EXISTENCE_AND_FULL_TORSION_CERTIFICATE_DONE' in s
assert not any(x in s for x in ['Runtime error','User error','syntax error','Assertion failed'])
for script,marker in [('tate_cusp_replay.py','C226_TATE_SYMMETRIC_PENCIL_PASS'),('groups_replay.py','C226_INDEPENDENT_TORSION_BOUND_PASS')]:
 with (out/(script+'.log')).open('w') as f:q=subprocess.run([sys.executable,str(P/script),str(out)],stdout=f,stderr=subprocess.STDOUT,cwd=ROOT)
 assert q.returncode==0 and marker in (out/(script+'.log')).read_text()
(out/'completion.json').write_text(json.dumps(dict(native_exit_code=r.returncode,complete=True,elapsed_seconds=time.monotonic()-start),indent=2)+'\n')
print('C226_OCCURRENCE_FRESH_REPRODUCTION_PASS '+str(out))
