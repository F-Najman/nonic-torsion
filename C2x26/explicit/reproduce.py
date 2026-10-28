from pathlib import Path
import subprocess,sys,json,time,argparse
P=Path(__file__).resolve().parent;ROOT=P.parents[1]
ap=argparse.ArgumentParser();ap.add_argument('output',type=Path);a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
assert not any(c in str(out) for c in ['"','\n','\r'])
wrapper=out/'driver.m';wrapper.write_text('C226Out:='+json.dumps(str(out))+';\nload '+json.dumps(str(P/'verify_native_reproduction.m'))+';\n')
t=time.monotonic()
with (out/'native.log').open('w') as f:r=subprocess.run(['magma','-b',str(wrapper)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
assert r.returncode==0 and 'C226_STANDALONE_EXPLICIT_OCCURRENCE_PASS' in (out/'native.log').read_text()
assert (out/'native_arithmetic.txt').read_bytes()==(ROOT/'C2x26/data/explicit/native_arithmetic.txt').read_bytes()
(out/'completion.json').write_text(json.dumps(dict(native_exit_code=r.returncode,complete=True,elapsed_seconds=time.monotonic()-t),indent=2)+'\n')
print('C226_EXPLICIT_FRESH_NATIVE_REPRODUCTION_PASS '+str(out))
