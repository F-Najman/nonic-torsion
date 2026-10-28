#!/usr/bin/env python3
"""Independent read-only Sage/M4RI ranks of retained binary matrices."""
from pathlib import Path
import gzip,json,time
from sage.all import GF,matrix
H=Path(__file__).resolve().parents[1]/'data'/'quotient45';t=time.monotonic();cases=json.loads(gzip.open(H/'canonical_ring_matrices.json.gz','rt').read());assert len(cases)==550
for c in cases:
 r=[int(a,16)for a in c['rows']];n=c['columns'];assert len(r)==132 and n>=232 and(n-56)%8==0;assert matrix(GF(2),[[a>>j&1 for j in range(n)]for a in r]).rank()==132
print('NONIC67H3_M4RI_READONLY_PASS cases550 rank132 elapsed',time.monotonic()-t)
