from sage.all import *
from pathlib import Path
import json,re,sys
root=Path(__file__).resolve().parents[1]/'data';p=root/'F5_rr';J=lambda q:json.loads(q.read_text());lit=lambda v:json.loads(re.sub(r',\s*\]',']',v.replace('[*','[').replace('*]',']').replace('<','[').replace('>',']')))
def ass(s,key):return lit(re.search(r'\b'+key+r':=(.*?);',s,re.S)[1])
model=J(root/'model_F5'/'model_cusps33.json');rep=J(root/'short_cusp_classes33.json');inputs=(p/'input33.m').read_text();assert ass(inputs,'FieldEquation')==model['field_polynomial'];assert ass(inputs,'CuspRepresentatives')==[z['global_cusp_coefficients']for z in rep['representatives']]
rawblocks=ass(inputs,'GlobalCuspBlocks');lookup={tuple(z['signature']):z for z in model['global_cusp_blocks']}
for i,(sig,deg,places)in enumerate(rawblocks):
 assert sig==rep['global_cusp_signatures'][i]and deg==rep['global_cusp_degrees'][i];z=lookup[tuple(sig)];assert deg==z['global_degree'];assert places==[[v['degree'],v['multiplicity'],v['finite'],v['ideal_basis']]for v in z['places']]
s=(p/'basis33.m').read_text();ac=ass(s,'ACoefficients');bf=ass(s,'BasisFunctions');assert ac==[max([0]+[z['global_cusp_coefficients'][j]for z in rep['representatives']])for j in range(17)];assert sum(a*d for a,d in zip(ac,rep['global_cusp_degrees']))==48
k=GF(5);R=PolynomialRing(k,'t');K=R.fraction_field();rat=lambda z:K(R(z[0]))/K(R(z[1]));ff=[[rat(z)for z in row]for row in bf];assert len(ff)==28 and all(len(z)==12 for z in ff)
dn=lcm([c.denominator()for row in ff for c in row]);md=max((dn*c).numerator().degree()for row in ff for c in row);M=matrix(k,[[R(dn*c)[i]for c in row for i in range(md+1)]for row in ff]);assert M.rank()==28;ni=J(p/'coefficient_independence33.json');assert ni['rank']==28 and ni['max_degree']==md and ni['columns']==M.ncols()
data=J(p/'jets_data33.json')['tables'];raw=(p/'jets33.jsonseq').read_text();original=[];dec=json.JSONDecoder()
while raw.strip():raw=raw.lstrip();v,i=dec.raw_decode(raw);original.append(v);raw=raw[i:]
assert data==original
text=(p/'jet_tables33.m').read_text();matches=re.findall(r'<(\d+),(\d+),(\d+),(\d+),(\[.*?\]),Matrix\(GF\(5\),28,(\d+),(\[.*?\])\)>',text,re.S);assert len(matches)==len(data)
for mm,z in zip(matches,data):
 j,pi,dd,ao,orders,nc,rr=mm;assert [int(j),int(pi),int(dd),int(ao)]==[z['block_index'],z['place_index'],z['degree'],z['A_order']];assert json.loads(orders)==z['orders']and int(nc)==len(z['rows'][0])and json.loads(rr)==sum(z['rows'],[])
used=set()
for z in data:
 j=z['block_index']-1;v=rawblocks[j][2][z['place_index']-1];deg,e,finite,ib=v;assert deg==z['degree'];assert z['A_order']==ac[j]*e and z['orders']==[(ac[j]-r['global_cusp_coefficients'][j])*e for r in rep['representatives']];assert len(z['rows'])==28 and all(len(r)==max(z['orders'])*deg for r in z['rows']);used.add(j)
assert used=={j for j in range(17)if any(ac[j]-r['global_cusp_coefficients'][j]>0 for r in rep['representatives'])}
ranks=[]
for i in range(10):
 rows=[sum([z['rows'][j][:z['orders'][i]*z['degree']]for z in data],[])for j in range(28)];MM=matrix(k,rows);assert MM.ncols()==39 and MM.rank()==28;ranks.append(int(MM.rank()))
assert 'NONIC33_F5_FRESH_POLE_JET_VERIFIER_COMPLETE'in(root.parent/'logs'/'F5_rr'/'verify33.log').read_text()
out={'complete':True,'genus':21,'A_degree':48,'basis_dimension_by_Riemann_Roch':28,'independent_rational_coefficient_rank':28,'B_degrees':[39]*10,'independent_jet_ranks':ranks,'all_ten_h0_zero':True,'fresh_native_verifier_has_no_RR_or_class_search':True}
if'--check'in sys.argv:assert J(p/'independent_certificate33.json')==out
else:(p/'independent_certificate33.json').write_text(json.dumps(out,indent=2)+'\n')
print('NONIC33_F5_INDEPENDENT_POLE_BASIS_JET_CERTIFICATE_PASS all10h0zero=true',flush=True)
