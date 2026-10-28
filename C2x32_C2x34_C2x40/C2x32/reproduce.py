from pathlib import Path
import subprocess,sys,json,time,resource,argparse,shutil
P=Path(__file__).resolve().parent;ROOT=P.parents[1]
ap=argparse.ArgumentParser();ap.add_argument('output',type=Path);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
assert not any(c in str(out) for c in ['"','\n','\r'])
wrapper=out/'driver.m';wrapper.write_text('C232Out:='+json.dumps(str(out))+';\nload '+json.dumps(str(P/'four_quotients_reproduction.m'))+';\n')
def limits():resource.setrlimit(resource.RLIMIT_CPU,(14400,14430))
start=time.monotonic()
with (out/'native.log').open('w') as f:r=subprocess.run(['magma','-b',str(wrapper)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits)
s=(out/'native.log').read_text();assert r.returncode==0 and 'C232_FOURQUOTIENT_DONE tested=107 remaining=0 complete=true' in s
assert not any(x in s for x in ['Runtime error','User error','syntax error','Assertion failed'])
assert (out/'four_quotients_arithmetic.txt').read_bytes()==(ROOT/'C2x32_C2x34_C2x40/data/C2x32/four_quotients_arithmetic.txt').read_bytes()
(out/'completion.json').write_text(json.dumps(dict(native_exit_code=r.returncode,complete=True,elapsed_seconds=time.monotonic()-start),indent=2)+'\n')
print('C232_EXCLUSION_FRESH_REPRODUCTION_PASS '+str(out))
