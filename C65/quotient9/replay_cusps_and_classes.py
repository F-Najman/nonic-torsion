from pathlib import Path
from math import gcd
from collections import Counter,deque
import ast,json,sys
P=Path(__file__).resolve().parents[1]/'data'/'quotient9';N=65
U=[u for u in range(N)if gcd(u,N)==1]
def chi(u):return(pow(u%5,2,5)in[1])==(pow(u%13,6,13)in[1])
H=[u for u in U if chi(u)];assert len(H)==24

def key(a,c):
 d=gcd(c,N);a%=d;return min((a,c%N),((-a)%d,(-c)%N))
C={key(a,c)for c in range(N)for a in range(gcd(c,N))if gcd(a,gcd(c,N))==1};assert len(C)==96

def dia(p,u):return key(p[0]*pow(u,-1,N),p[1]*u)
def orbH(p):return frozenset(dia(p,u)for u in H)
YC=set(map(orbH,C));assert len(YC)==8;index={p:y for y in YC for p in y}
def gal(p,u):return key(p[0]*u,p[1])
def yg(y,u):return index[gal(min(y),u)]
def orbits(C,act,us):
 todo=set(C);ans=[]
 while todo:
  p=min(todo,key=repr);o=frozenset(act(p,u)for u in us);assert o<=todo;todo-=o;ans.append(o)
 return ans
F=[];v=1
while v not in F:F.append(v);v=v*2%65
SG=orbits(C,gal,U);SF=orbits(C,gal,F);YG=orbits(YC,yg,U);YF=orbits(YC,yg,F)
assert Counter(map(len,SG))=={1:24,4:6,12:2,24:1}and Counter(map(len,SF))=={1:24,4:6,6:4,12:2}
small=[o for o in SF if len(o)<=3];assert all(o in SG for o in small)and len(small)==24
bad=[o for o in SF if len(o)<=9 and o not in SG];assert len(bad)==4 and all(len(o)==6 for o in bad)
for o in bad:
 ys={index[p]for p in o};assert len(ys)==1 and N//gcd(min(o)[1],N)==1;yy=next(iter(ys));assert frozenset([yy])in YF
assert Counter(map(len,YG))=={1:2,2:3}and Counter(map(len,YF))=={1:4,2:2}
# The dual j-width is N/(point j-width); global Galois action is unchanged.
qrows=[]
for w in[1,5,13,65]:
 wg=[o for o in YG if gcd(min(next(iter(o)))[1],N)==w];wf=[o for o in YF if gcd(min(next(iter(o)))[1],N)==w];qrows.append(dict(dual_width=w,global_degrees=sorted(map(len,wg)),finite_degrees=sorted(map(len,wf))))
assert qrows==[dict(dual_width=1,global_degrees=[1,1],finite_degrees=[1,1]),dict(dual_width=5,global_degrees=[2],finite_degrees=[2]),dict(dual_width=13,global_degrees=[2],finite_degrees=[2]),dict(dual_width=65,global_degrees=[2],finite_degrees=[1,1])]

def rd(f):return ast.literal_eval((P/f).read_text().replace('[*','[').replace('*]',']').replace('<','[').replace('>',']').replace('false','False').replace('true','True'))
g,rels,cd,widths,sig,sf,hc,ac,jci,ai,kg,aop=rd('finite_group_A3_data.txt');base,tests=rd('non_global_cusp_tests.txt');assert g==[4,16,336,0]and ai==[4,8]and base==1 and widths[base-1]==1
mods=g[:-1]
def ky(v):assert v[-1]==0;return tuple(a%b for a,b in zip(v[:-1],mods))
def add(a,b):return tuple((v+w)%m for v,w,m in zip(a,b,mods))
def mul(n,a):return tuple(n*v%m for v,m in zip(a,mods))
hc=list(map(ky,hc));ac=list(map(ky,ac));aop=list(map(ky,aop));zero=(0,0,0)
assert cd==[1,1,1,2,2,1]and sig==[6,3,2,4,5,1]and sf==[1,3,2,4,5,6]
assert hc[base-1]==zero
for i in range(6):
 expected=add(add(hc[sig[i]-1],mul(-1,hc[i])),mul(-cd[i],add(hc[sig[base-1]-1],mul(-1,hc[base-1]))))
 assert ac[i]==expected
# Reconstruct the dual A3 boundary coefficients and every base correction.
def rawA(i):
 row=[0]*6
 for j,n in[(sig[sf[i]-1]-1,1),(sf[i]-1,3),(sig[i]-1,-1),(i,-3)]:row[j]+=n
 return row
rows=[rawA(i)for i in range(6)];baseA=rows[base-1]
for i in range(6):
 expected=zero
 for j in range(6):expected=add(expected,mul(rows[i][j]-cd[i]*baseA[j],ac[j]))
 assert expected==aop[i]
operator={zero:zero};todo=deque([zero])
while todo:
 v=todo.popleft()
 for a,b in zip(ac,aop):
  vv=add(v,a);ww=add(operator[v],b)
  if vv in operator:assert operator[vv]==ww
  else:operator[vv]=ww;todo.append(vv)
assert len(operator)==32;kernel={v for v,w in operator.items()if w==zero};assert len(kernel)==8
assert all(sig[sig[i]-1]==i+1 for i in range(6))
assert all(sf[i]==(sig[i]if widths[i]==65 else i+1)for i in range(6))
expected=[i+1 for i,w in enumerate(widths)if w==65];assert[len(v)for v in tests]==[2,2]and[i for i,b in tests]==expected
for i,b in tests:assert not b and mul(6,ac[i-1])not in kernel
out=dict(complete=True,level=65,prime=2,source_geometric_cusps=96,source_global_signature=dict(Counter(map(len,SG))),source_F2_signature=dict(Counter(map(len,SF))),source_small_cusps_globally_intact=True,non_global_sextic_cusps=4,all_project_to_dual_width65_with_norm6=True,quotient_genus=9,quotient_geometric_cusps=8,quotient_width_rows=qrows,finite_J_invariants=g[:-1],finite_alpha_cusp_image_order=32,mu_A3_boundary_and_base_corrections_independently_recomputed=True,mu_A3_kernel_order=8,both_split_cusp_norm6_classes_excluded=True,all_non_global_allbad_reductions_excluded=True)
if '--readonly'in sys.argv:assert json.loads((P/'cusp_class_sage_certificate.json').read_text())==json.loads(json.dumps(out))
else:(P/'cusp_class_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C65_NONGLOBAL_CUSP_CLASS_INDEPENDENT_PASS',json.dumps(out),flush=True)
