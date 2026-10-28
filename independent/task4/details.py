# Per-curve details over F_{2^e} for selected N, from the detail files written by hist.gp:
# number of F_q-isomorphism classes E with N | exponent, their group orders, j-degrees over F_2, and j=0 (supersingular) cases.
import glob, re, os, collections
D = os.path.dirname(os.path.abspath(__file__))
pat = re.compile(r'\[(\d+), (\d+), "(\w+)", (\d+), \[(\d+), (\d+)\], ([\d/]+), (\d+), "(.*)"\]')
rows = collections.defaultdict(list)
for fn in glob.glob(os.path.join(D, 'hist', 'd_2_*.txt')):
    for line in open(fn):
        m = pat.match(line.strip())
        if not m: continue
        p, e, tag, N, n1, n2, w, jdeg, j = m.groups()
        rows[(int(N), int(e))].append((tag, int(n1)*int(n2), (int(n1), int(n2)), w, int(jdeg)))
for N in [27, 39, 53, 67, 69, 71, 73, 75, 79, 83, 85, 87, 89, 91, 93]:
    for e in range(1, 10):
        r = rows.get((N, e), [])
        if not r: continue
        c = collections.Counter((tag, order, struct, jdeg, w) for tag, order, struct, w, jdeg in r)
        print('N=%d F_2^%d: %d classes; ' % (N, e, len(r)) + '; '.join('%dx[%s #E=%d %s jdeg=%d w=%s]' % (v, k[0], k[1], k[2], k[3], k[4]) for k, v in sorted(c.items())))
