#!/usr/bin/env python3
"""Read-only independent M4RI rank replay of explicit small certificates."""
from pathlib import Path
import gzip,json,time
from sage.all import GF,matrix
H=Path(__file__).resolve().parents[1]/'data'/'quotient19';t0=time.monotonic()
with gzip.open(H/'canonical_ring_matrices.json.gz','rt')as fp:cases=json.load(fp)
assert len(cases)==32
F=GF(2)
for c in cases:
 rows=[int(x,16)for x in c['rows']];assert len(rows)==90 and c['columns']==164
 A=matrix(F,[[x>>j&1 for j in range(164)]for x in rows]);assert A.rank()==90
print('NONIC73_SMALL_MATRICES_SAGE_M4RI_PASS cases=32 rank=90 time=',time.monotonic()-t0)
