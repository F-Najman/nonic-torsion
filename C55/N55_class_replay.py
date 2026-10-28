from pathlib import Path
import ast
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
H=Path(__file__).resolve().parent
s=(H/'data'/'N55_quotient_cusp_arithmetic.txt').read_text().replace('[*','[').replace('*]',']')
inv,rels,gens,hinv,degrees,widths=ast.literal_eval(s)
assert inv==[5,5,10,1200,0] and hinv==[5,10,1200]
assert sorted(degrees)==[1]*4+[2]*2+[4]*2 and all(r[-1]==0 for r in rels+gens)
relations=[r[:-1]for r in rels if any(r)];assert len(relations)==4
assert abs(Matrix(relations).det())==300000
lattice=hermite_normal_form(Matrix(relations+[r[:-1]for r in gens]).T)
index=abs(int(lattice.det()));order=300000//index
assert index==5 and order==60000 and order%32==0 and order%64!=0
print('NONIC55_CUSPIDAL_IMAGE_HNF_REPLAY_PASS finite_J_order=300000 rational_cuspidal_image_order=60000 two_primary_order=32')
