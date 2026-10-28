"""Read-only check that the fresh degree9 native run finished all lower branches."""
from pathlib import Path
import re,ast,sys
P=Path(__file__).resolve().parents[1]
s=Path(sys.argv[1] if len(sys.argv)>1 else P/'logs/N26_lower/fresh_D9_completed_lower8_prefix.log').read_text()
a=ast.literal_eval((P/'data/degree8/N26_D8_class_arithmetic.txt').read_text().replace('[*','[').replace('*]',']').replace('false','False').replace('true','True'))
inv,rel,gens,cubic,records,fp,nc,reps=a
m=re.search(r'NONIC_N26_CLASS_GROUP (\[[^\]]*\])\s*(\[[^\]]*\])\s*(\[\s*(?:\[[^\]]*\],?\s*){6}\])',s)
assert m and ast.literal_eval(m.group(1))==inv and ast.literal_eval(m.group(2))==[2,14,266,3990,11970] and ast.literal_eval(m.group(3))==cubic
rs=re.findall(r'NONIC_N26_CLASS_TEST ([5-8]) (\d+) false \[\]\s*(\[[^\]]*\])',s)
assert len(rs)==239
for (ds,ii,vs),old in zip(rs,records):
    assert [int(ds),int(ii)]==old[:2] and ast.literal_eval(vs)==old[2] and old[3] is False and old[4]==[]
for d,n in zip(range(5,9),[6,26,48,159]):assert f'NONIC_N26_CLASS_TYPE {d} {n} 0 ' in s
assert not any(x in s for x in ['Runtime error','User error','syntax error','Assertion failed'])
print('N26_D9_FRESH_COMPLETED_LOWER8_PREFIX_PASS representatives=239 cubic_shift_degrees=5,6')
