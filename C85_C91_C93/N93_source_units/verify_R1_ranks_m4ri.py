#!/usr/bin/env python3
"""Independent Sage/M4RI ranks, without native launch or output mutation."""
from pathlib import Path
import gzip,json,time
import os
os.environ["OPENBLAS_NUM_THREADS"]="1";os.environ["OMP_NUM_THREADS"]="1"
from sage.all import GF,matrix
H=Path(__file__).resolve().parents[1]/'data/N93_source_units';t=time.monotonic();rr=json.loads(gzip.open(H/'canonical_R1_matrices.json.gz','rt').read());assert len(rr)==120
for c in rr:
 rows=[int(v,16)for v in c['rows']];n=c['columns'];assert len(rows)==261 and n>=519 and all(v.bit_length()<=n for v in rows);assert matrix(GF(2),[[v>>j&1 for j in range(n)]for v in rows]).rank()==261
print('NONIC93_R1_M4RI_RANK_PASS cases120 rank261 elapsed',time.monotonic()-t)
