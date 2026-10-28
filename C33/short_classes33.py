from sage.all import *
from pathlib import Path
import json,re,collections
p=Path(__file__).resolve().parent/'data';d=p/'rebound';binding=json.loads((d/'cusp_binding33.json').read_text());source=(d/'class_arithmetic33.py').read_text()
def assign(name):
 m=re.search(r'\b'+name+r'\s*=\s*(.*?)(?=\n\w+\s*=|\Z)',source,re.S);assert m,name
 return json.loads(re.sub(r',\s*\]',']',m[1].strip().replace('True','true').replace('False','false')))
qc=assign('rational_cusp_generators');assert qc==binding['global_cusp_classes'];mods=binding['Pic_invariants'][:2];base=binding['base_cusp_index']-1;assert qc[base]==[0,0,0]
known=[assign(z)for z in ['known_septic_plus_degree2cusp','known_octic_plus_cusp','known_nonic_classes']];knownset={tuple(z)for vv in known for z in vv};qe={tuple(z)for z in assign('rational_effective_cusps')};tests=assign('tests');candidates={tuple(v)for deg,v,h,q in tests if h and not q};unknown=sorted(candidates-knownset);assert len(unknown)==10 and all(next(deg for deg,v,h,q in tests if tuple(v)==z)==9 for z in unknown)
idx=[i for i in range(17)if i!=base];deg=binding['global_cusp_degrees'];A=matrix(ZZ,[qc[i][:2]for i in idx]+[[mods[0],0],[0,mods[1]]]);D,U,V=A.smith_form();assert U*A*V==D;ker=A.transpose().right_kernel().basis_matrix()[:,:16];assert ker.rank()==16;weights=diagonal_matrix(ZZ,[deg[i]for i in idx]);B=(ker*weights).LLL();GS,mu=B.gram_schmidt();reps=[]
def cuspc(w):
 co=[ZZ(0)]*17
 for i,a in zip(idx,w):co[i]=ZZ(a)
 co[base]=9-sum(co[i]*deg[i]for i in idx);assert sum(co[i]*deg[i]for i in range(17))==9;return co
def obj(w):
 co=cuspc(w);return(sum(max(a,0)*di for a,di in zip(co,deg)),sum(abs(a)*di for a,di in zip(co,deg)))
for target in unknown:
 tv=vector(ZZ,target[:2])*V;xx=vector(ZZ,18)
 for j in range(2):assert tv[j]%D[j,j]==0;xx[j]=tv[j]//D[j,j]
 particular=vector(ZZ,(xx*U)[:16]);rest=vector(QQ,particular*weights)
 for j in range(15,-1,-1):n=ZZ((rest.dot_product(GS[j])/GS[j].dot_product(GS[j])).round());rest-=n*B[j]
 w=vector(ZZ,[rest[j]/deg[idx[j]]for j in range(16)]);kb=[vector(ZZ,[row[j]/deg[idx[j]]for j in range(16)])for row in B.rows()]
 for _ in range(200):
  best=min([w]+[w+sgn*r for r in kb for sgn in [-1,1]],key=obj)
  if obj(best)>=obj(w):break
  w=best
 co=cuspc(w);got=[sum(co[i]*qc[i][j]for i in range(17))%mods[j]for j in range(2)];assert got==list(target[:2])
 reps.append({'finite_class':list(target),'global_cusp_coefficients':list(map(int,co)),'degree':9,'positive_degree':int(obj(w)[0]),'negative_degree':int(obj(w)[0]-9)})
out={'complete':True,'base_cusp_index':base+1,'base_cusp_signature':binding['base_signature'],'global_cusp_degrees':deg,'global_cusp_signatures':binding['global_cusp_signatures'],'Pic_invariants':binding['Pic_invariants'],'unknown_classes':10,'representatives':reps};(p/'short_cusp_classes33.json').write_text(json.dumps(out,indent=2)+'\n');print('NONIC33_SHORT_GLOBAL_CUSP_REPRESENTATIVES_COMPLETE',[(z['positive_degree'],z['negative_degree'])for z in reps],flush=True)
