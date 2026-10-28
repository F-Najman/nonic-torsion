"""Recompute the complete nonic C2x30 exclusion in a NEW directory.

Run from any directory: python3 reproduce_D9.py NEW_RUN_DIRECTORY.
Requires Magma and Python with SymPy. Every native stage has a600s limit;
no previous-degree arithmetic result or retained session is loaded.
"""
from pathlib import Path
import argparse,subprocess,sys,json,shutil,time

P=Path(__file__).resolve().parent
ROOT=P.parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);ap.add_argument('--magma',default='magma');args=ap.parse_args()
    out=args.run.resolve();out.mkdir(parents=True,exist_ok=True)
    assert not any(c in str(out) for c in ['"','\n','\r'])
    shutil.copy2(P/'logs/run/exceptional_pencil_inputs.m',out/'exceptional_pencil_inputs.m')
    stages=[('classes256_D9','A3_D9_RR_PROBE_DONE'),('exceptional_fibres_D9','EXCEPTIONAL_FIBRES_DONE'),('three_quotients_D9','THREE_QUOTIENTS_DONE')]
    records=[]
    for stage,marker in stages:
        wrapper=out/(stage+'_driver.m');pre=''
        if stage=='exceptional_fibres_D9':pre='load '+json.dumps(str(out/'exceptional_pencil_inputs.m'))+';\n'
        if stage=='three_quotients_D9':pre='load '+json.dumps(str(out/'exceptional_fibres_data.m'))+';\n'
        wrapper.write_text('C230Out:='+json.dumps(str(out))+';\n'+pre+'load '+json.dumps(str(P/(stage+'.m')))+';\n')
        console=out/(stage+'_console.log');start=time.monotonic()
        with console.open('w') as f:r=subprocess.run([args.magma,'-b',str(wrapper)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,timeout=600)
        txt=console.read_text();assert r.returncode==0 and marker in txt
        assert not any(z in txt for z in ['Runtime error','User error','syntax error','Assertion failed'])
        records.append(dict(stage=stage,returncode=r.returncode,elapsed_seconds=time.monotonic()-start))
        print('C230_D9_FRESH_NATIVE_PASS '+stage,flush=True)
    with (out/'independent_replay.log').open('w') as f:r=subprocess.run([sys.executable,str(P/'replay_D9.py'),str(out)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,timeout=600)
    assert r.returncode==0 and 'C230_D9_INDEPENDENT_REPLAY_PASS' in (out/'independent_replay.log').read_text()
    (out/'completion.json').write_text(json.dumps(dict(degree=9,complete=True,excluded_groups=[[2,30],[1,60]],native_stages=records),indent=2)+'\n')
    print('C230_D9_FRESH_REPRODUCTION_PASS '+str(out))

if __name__=='__main__':main()
