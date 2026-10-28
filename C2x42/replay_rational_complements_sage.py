from pathlib import Path
import sys,runpy,json,contextlib,io
CODE=Path(__file__).resolve().parent;P=CODE/'data';readonly='--readonly'in sys.argv;argv=sys.argv[:]
sys.argv=[str(CODE/'genus5_quotient/replay_genus5_sage.py'),'--readonly']
with contextlib.redirect_stdout(io.StringIO()):g=runpy.run_path(str(CODE/'genus5_quotient/replay_genus5_sage.py'))
sys.argv=argv
K=g['KQ'];hc=g['hc'];tor=g['branches'][0];nr=K.key([3*v for v in hc[tor[0]]]);zero=(0,0,0,0);assert nr!=zero and all(K.key([3*v for v in hc[j]])==nr for j in tor)and all(K.key([6*v for v in hc[j]])==zero for j in tor)
a=g['rd'](P/'genus5_quotient/filter_919/filter_arithmetic.txt');cases=[[1,4,5],[1,6,39],[2,4,3,4,55],[2,4,3,4,65],[2,4,3,5,260],[2,4,3,5,461],[2,4,3,5,589]];outrows=[]
for d in cases:
 r=[r for r in a[8]if r[1]==d];assert len(r)==1;r=r[0];assert r[3]and r[5]and K.key(r[4])==zero and 0<=r[2]<=5;outrows.append(dict(descriptor=d,complement_degree=r[2],normalized_image_class=r[4],coset=list(zero)))
cusp=json.loads((P/'CM_Y5_retest/source_cusp_certificate.json').read_text());assert cusp['global_signature']=={'1':18,'2':9,'6':6}and cusp['F5_signature']=={'1':18,'2':9,'3':6,'6':3}and cusp['all_degree_le_two_cusps_globally_intact']
cm=g['rd'](P/'CM_Y5_retest/retests.txt');assert cm[1]==[1,2,6,8]and all(r[4]==[0]for ii in cm[1]for r in cm[0][ii-1][3]);assert json.loads((P/'CM_Y5_retest/sage_certificate.json').read_text())['complete']
out=dict(complete=True,source_cusps_below_degree6=[1,2,3],only_nonrational_option_single_split_cubic=True,split_cubic_image_coset=list(nr),split_cubic_image_order=2,nonCM_cases=outrows,CM_cases=[1,2,6,8],all_11_complements_globally_rational=True)
if readonly:assert json.loads((P/'rational_complement_certificate.json').read_text())==out
else:(P/'rational_complement_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C242_ALL_FINAL_CUSP_COMPLEMENTS_GLOBALLY_RATIONAL_SAGE_PASS',json.dumps(out),flush=True)
