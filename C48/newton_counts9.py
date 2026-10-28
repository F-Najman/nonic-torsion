#!/usr/bin/env python3
"""Numbers of closed points of degree 1..9 of X_1(48) (source) and Y = X_1(2,24) (target) over F_5, from the
Frobenius polynomials in the degree-8 target setup export (Newton identities); writes C48/data/degree9_Newton_counts.json.
Extracted for this repository from the inline Python command of the Codex session (transcript 2026-10-26T19:40:08Z)
that wrote degree9_Newton_counts.json; only the paths are adapted (the clean_run_v2 copy of
target_2_24_arithmetic.py is byte-identical to the clean_run copy read here)."""
from pathlib import Path
import ast,json
p=Path(__file__).resolve().parent;o=p/'data'/'degree8';D={n.targets[0].id:ast.literal_eval(n.value) for n in ast.parse((o/'clean_run/target_2_24_arithmetic.py').read_text()).body}
def closed(coeff,p,n):
 a=list(reversed(coeff));ss=[0];out=[0]
 for m in range(1,n+1):
  sm=-m*a[m]-sum(a[j]*ss[m-j] for j in range(1,m));ss.append(sm);nm=p**m+1-sm;cm=nm-sum(d*out[d] for d in range(1,m) if m%d==0);assert cm%m==0;out.append(cm//m)
 return out
source=closed(D['source_frobenius_coefficients'],5,9);target=closed(D['target_frobenius_coefficients'],5,9);print(source,target)
(p/'data'/'degree9_Newton_counts.json').write_text(json.dumps({'source_closed1to9':source[1:],'target_closed1to9':target[1:],'source_frobenius_coefficients':D['source_frobenius_coefficients'],'target_frobenius_coefficients':D['target_frobenius_coefficients']},indent=2)+'\n')
