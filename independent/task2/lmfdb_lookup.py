# polredabs each realization field (from the entry logs) and look it up in the LMFDB (devmirror).
import sys, re, glob, os
sys.path.append('/home/fnajman/.local/lib/python3.11/site-packages')
import cypari2, psycopg2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
D = os.path.dirname(os.path.abspath(__file__))
con = psycopg2.connect(host='devmirror.lmfdb.xyz', dbname='lmfdb', user='lmfdb', password='lmfdb', connect_timeout=30)
cur = con.cursor()
for fn in sorted(glob.glob(os.path.join(D, 'entry_*.log'))):
    txt = open(fn).read()
    m = re.search(r'FIELDPOLY (.*)', txt)
    head = re.search(r'ENTRY (.*)', txt).group(1)
    f = pari(m.group(1).replace('t', 'x'))
    g = pari.polredabs(f)
    nf = pari.nfinit(g)
    disc = int(pari.nfdisc(g))
    coeffs = [int(c) for c in pari.Vecrev(g)]
    cur.execute("select label, coeffs, disc_abs, disc_sign, galois_label, class_number from nf_fields where degree=9 and disc_abs=%s", (abs(disc),))
    rows = cur.fetchall()
    hit = None
    for r in rows:
        if [int(c) for c in r[1]] == coeffs: hit = r
    iso = []
    if hit is None:
        for r in rows:
            h = pari('x^0') * 0
            gg = pari('Pol(Vecrev(%s))' % [int(c) for c in r[1]])
            if pari.nfisisom(g, gg): iso.append(r[0])
    print(os.path.basename(fn), '|', head)
    print('   polredabs', g, ' disc', disc, ' LMFDB rows with this |disc|:', len(rows))
    print('   LMFDB label:', hit[0] if hit else ('isomorphic: %s' % iso if iso else 'NOT IN LMFDB'), ('galois %s class_number %s' % (hit[4], hit[5])) if hit else '')
