#!/usr/bin/env python3
"""Independent exact arithmetic for the nonic prime bound and binary FI matrices.
No Magma packages are imported. The cited mathematical theorems and the
native construction of the exported Hecke matrices are separate proof inputs.
"""
from pathlib import Path
from fractions import Fraction
import json,re,hashlib,math
from sympy import primerange,n_order,primitive_root
D=Path(__file__).resolve().parent

def rank(rows):
    piv={}
    for row in rows:
        while row:
            p=row.bit_length()-1
            if p in piv:row^=piv[p]
            else:piv[p]=row;break
    return len(piv)

def compositions(parts,total):
    if parts==1:
        yield (total,);return
    for a in range(total+1):
        for b in compositions(parts-1,total-a):yield(a,)+b

# LaTeX source (alttext) of Table 2 of the HTML version of Derickx-Stoll, arXiv:2505.14109v1
t=(D/'data'/'derickx_stoll_table2.txt').read_text()
ps=re.findall(r'(?<![A-Za-z])p&([^\\]+)',t)
orders=re.findall(r'\\ord\(\\chi\)&([^\\]+)',t)
assert len(ps)==len(orders)==9
strange={}
for pp,oo in zip(ps,orders):
    pp=[int(x)for x in pp.split('&') if x.strip()]
    oo=[math.lcm(*map(int,x.split(','))) for x in oo.split('&') if x.strip()]
    assert len(pp)==len(oo)
    strange.update(zip(pp,oo))
assert len(strange)==74
# Compare independently with the native authors' table snapshot.
a=(D/'strange_primes_original.magma').read_text().split('TABLE2 := [',1)[1].split('];',1)[0]
orig={}
for p,o,dim in re.findall(r'<(\d+),\s*(\d+),\s*(\d+)>',a):orig[int(p)]=math.lcm(orig.get(int(p),1),int(o))
assert strange==orig
# (3^(9/2)+1)^2=19684+162 sqrt(3) lies strictly between19964 and19965.
assert 280**2 < 162**2*3 <281**2
cap=19964
assert 325*(127**2-1)>65536*8*9
excluded=[]
large_data=[]
for p in primerange(127,cap+1):
    st=strange.get(p,1)
    assert (p-1)%(2*st)==0
    h=(p-1)//(2*st)
    assert 9%h!=0
    large_data.append({'p':p,'strangeness':st,'H_order':h})
    excluded.append(p)
small_data=[]
for p,a,n in [(103,2,5),(107,2,5),(109,int(primitive_root(109)),6),(113,3,7)]:
    assert p not in strange
    o=int(n_order(a,p));o=o//2 if o%2==0 else o
    assert o==(p-1)//2
    g=Fraction(325*(p*p-1),65536)
    assert g>n*9
    assert 9%o!=0
    small_data.append({'p':p,'diamond_generator':a,'diamond_order':o,'n':n,'gonality_lower_bound':[g.numerator,g.denominator]})
    excluded.append(p)
fi=[]
for p in [97,101]:
    lines=(D/'data'/f'FI_prime{p}_rows.txt').read_text().splitlines()
    pp,ind,n,r,nrows,ncols=map(int,lines[0].split())
    assert pp==p and len(lines)==nrows+1 and nrows==9+4*(ind-1)
    assert all(len(x)==ncols and set(x)<={'0','1'}for x in lines[1:])
    rows=[int(x,2)for x in lines[1:]]
    count=0
    for sizes in compositions(ind,9):
        if sizes[0]!=max(sizes):continue
        chosen=rows[:sizes[0]]
        for k in range(1,ind):chosen+=rows[9+(k-1)*4:9+(k-1)*4+sizes[k]]
        assert len(chosen)==9 and rank(chosen)==9,(p,sizes)
        count+=1
    o=int(n_order(2,p));o=o//2 if o%2==0 else o
    assert o>9 # nonrational cusps cannot occur in the reduction.
    line=next(x for x in (D/'logs'/'finite_field_prime_orders.log').read_text().splitlines()if x.startswith(f'TORSION n={p} p=2 '))
    assert 'over F_p^k = [] ;'in line
    fi.append({'p':p,'index_H':ind,'Tn':n,'q':3,'diamond_generator':r,'rows':nrows,'columns':ncols,'binary_rank':rank(rows),'prefix_configurations':count,'all_ranks':9,'nonrational_cusp_reduction_degree':o})
    excluded.append(p)
possible=list(primerange(2,90))
assert sorted(set(possible)|set(excluded))==list(primerange(2,cap+1))
assert len(excluded)==len(set(excluded))
assert 'NONIC_PRIME_FI_CERTIFICATE_COMPLETE'in(D/'logs'/'prime_FI_certificate.log').read_text()
assert 'NONIC_PRIME_FINITE_FIELD_ENUMERATION_COMPLETE'in(D/'logs'/'finite_field_prime_orders.log').read_text()
res={'degree':9,'status':'PROVED PRIME UPPER BOUND; NOT FULL CLASSIFICATION','prime_bound':89,'possible_primes':possible,'known_occurring_primes':list(primerange(2,32)),'unresolved_prime_orders':[p for p in possible if p>31],'oesterle_integer_cap':cap,'all_primes_in_oesterle_range':len(possible)+len(excluded),'excluded_prime_orders':sorted(excluded),'large_prime_details':large_data,'small_global_exclusions':small_data,'formal_immersion_exclusions':fi,'strange_primes_in_oesterle_range':{str(p):v for p,v in strange.items()if p<=cap}}
(D/'data'/'prime_bound89_verified.json').write_text(json.dumps(res,indent=2)+'\n')
print('NONIC_PRIME_BOUND_PASS bound=89 possible='+str(len(possible))+' excluded='+str(len(excluded))+' unresolved_prime_orders='+str(res['unresolved_prime_orders']))
for x in fi:print('BINARY_FI_INDEPENDENT_REPLAY',x)
