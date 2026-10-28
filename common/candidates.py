#!/usr/bin/env python3
"""Audited minimal-obstruction bookkeeping for Phi(9), never a search-based classification."""
from pathlib import Path
import hashlib,json
from sympy import divisors,primefactors,primerange
H=Path(__file__).resolve().parent/'data';ROOT=H.parents[1]

def subgroups(g):
 return {(int(a),int(b)) for a in divisors(g[0]) for b in divisors(g[1]) if b%a==0}
def subgroup(a,b):return b[0]%a[0]==0 and b[1]%a[1]==0

def read_ledger():
 # The groups of the accepted proof ledger of the original record, with the folder of this
 # repository that holds each proof (common/data/ledger.json); the SHA-256 of the evidence
 # files are checked by the tests (common/data/evidence_sha256.tsv).
 d=json.loads((H/'ledger.json').read_text());assert d['degree']==9
 for k in ('realizations','exclusions'):
  for r in d[k]:
   assert r['degree']==9 and r['audited'] is True
   assert r['folder'] is None or (ROOT/r['folder']).is_dir(),r['folder']
 assert d['prime_bound']['audited'] is True
 return d

def frontier(base,primes):
 assert all(subgroups(g)<=base for g in base)
 cyclic={n for m,n in base if m==1}
 # A minimal missing cyclic n is prime, or has n/p in cyclic for every p|n.
 # Generating n=c*p, c in cyclic, is exhaustive and avoids a long interval.
 trial=set(primes)|{c*p for c in cyclic for p in primes}
 answer={(1,n) for n in trial if (1,n) not in base and
         all((1,n//p) in base for p in primefactors(n))}
 # For a minimal noncyclic obstruction its cyclic subgroup C_n is proper.
 for n in cyclic:
  g=(2,n)
  if n%2==0 and g not in base and subgroups(g)-{g}<=base:answer.add(g)
 return sorted(answer)

def main():
 d=read_ledger();primes=set(map(int,primerange(2,d['prime_bound']['upper_bound']+1)))
 inf={tuple(g) for g in d['published_infinite_families']['groups']}
 known=inf|{tuple(r['group']) for r in d['realizations']}
 excluded={tuple(r['group']) for r in d['exclusions']}
 assert len(inf)==42 and all(m in (1,2) for m,n in known)
 assert not any(subgroup(e,g) for e in excluded for g in known)
 f=frontier(known,primes);eliminated=[g for g in f if any(subgroup(e,g) for e in excluded)]
 remaining=[g for g in f if g not in eliminated]
 out={'degree':9,'status':'COMPLETE' if not remaining else 'INCOMPLETE RESEARCH',
 'classification_complete':not remaining,'prime_bound':d['prime_bound']['upper_bound'],
 'prime_universe':sorted(primes),'first_invariant_factors':[1,2],
 'phi_infinite':sorted(inf),'known_lower_bound':sorted(known),
 'classification':sorted(known) if not remaining else None,
 'minimal_outside_known':f,'certified_exclusions':sorted(excluded),
 'eliminated_minimal_obstructions':eliminated,'remaining_minimal_obstructions':remaining,
 'warning':('All minimal obstructions are excluded; the displayed51 groups give the complete classification.' if not remaining else 'Remaining minimal obstructions are proof tasks, not a list of occurring or possible full groups; an occurrence can expose larger obstructions.')}
 (H/'candidates.json').write_text(json.dumps(out,indent=2)+'\n')
 def label(g):return f'C{g[1]}' if g[0]==1 else f'C{g[0]} x C{g[1]}'
 lines=[f'Phi(9): {out["status"]}',f'Known occurring groups: {len(known)}',
 f'Prime universe: primes <= {out["prime_bound"]}',f'Minimal obstructions: {len(f)}',
 f'Excluded minimal obstructions: {len(eliminated)}',f'Remaining minimal obstructions: {len(remaining)}','',
 'Known occurring groups:',', '.join(map(label,sorted(known))),'',
 'Remaining minimal obstructions:',(', '.join(map(label,remaining)) or 'None'),'',out['warning'],'NONIC_CANDIDATE_AUDIT_PASS']
 (H/'candidates.txt').write_text('\n'.join(lines)+'\n');print('\n'.join(lines))
if __name__=='__main__':main()
