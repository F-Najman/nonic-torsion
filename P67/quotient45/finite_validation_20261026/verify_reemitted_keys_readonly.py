#!/usr/bin/env python3
"""Independent retained/fresh key binding and exact Frobenius fiber counts."""
from pathlib import Path
import ast,re,json,sys
H=Path(__file__).resolve().parents[2]/'data'/'quotient45'/'finite_validation_20261026';Q=H.parent

def lit(s,n):return ast.literal_eval(re.search(re.escape(n)+r':=(.*?);',s,re.S)[1].replace('[*','[').replace('*]',']').replace('<','(').replace('>',')'))
s=(H/'validated_keys.m').read_text();old=(Q/'source_places_data.m').read_text();assert lit(s,'ValidatedSourceKeys')==lit(old,'SourcePlaceKeys');assert lit(s,'ValidatedSourcePermutation')==lit(old,'SourceSigmaPermutation');assert lit(s,'ValidatedSourceOrbits')==lit(old,'SourceSigmaOrbits')
# Target lists are appended sequentially; parse each actual append.
kt=[ast.literal_eval(v.replace('[*','[').replace('*]',']').replace('<','(').replace('>',')'))for v in re.findall(r'Append\(~ValidatedTargetKeys,(.*?)[)]\s*;',s,re.S)]
pt=[ast.literal_eval(v)for v in re.findall(r'Append\(~ValidatedTargetPermutations,(.*?)[)]\s*;',s,re.S)];assert len(kt)==len(pt)==7
for i,(ks,pp)in enumerate(zip(kt,pt),1):
 t=(Q/f'target_places{i}_jbase_data.m').read_text();assert ks==lit(t,'TargetPlaceKeys')and pp==lit(t,'TargetSigmaPermutation')and len(ks)==11
G={1,66,38,29,37,30};assert len(G)==6 and {a*b%67 for a in G for b in G}==G and(67-1)//len(G)==11
counts=lit(s,'FrobeniusFiberCounts');assert len(counts)==12
for f,q,n in counts:assert q==pow(2,f,67)and q not in G and n==11
assert sorted({a[0]for a in counts})==[7,8,9,24,28,36]
print('NONIC67_FRESH_KEY_BINDING_READONLY_PASS source55 target77 fibers12 exact_count11 all_raw_and_reemitted_keys_identical=true')
