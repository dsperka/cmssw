#!/usr/bin/env python3
"""What the legacy refit (with the MTD hit) does to the momentum, and legacy vs Alpaka inputs of the time of flight.

usage: python3 analyze_refit.py table.npz      (table from make_table.py)
"""
import sys

import numpy as np

z = np.load(sys.argv[1] if len(sys.argv) > 1 else 'table.npz')
c = {n: i for i, n in enumerate(z['cols'])}
d = z['data']
col = lambda n: d[:, c[n]]
p, sp, pref, spref, ptrue = col('p'), col('sp'), col('pref'), col('spref'), col('ptrue')
region = np.where(np.abs(col('eta')) < 1.5, 'BTL', 'ETL')
bins = [(0.7, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 5.0), (5.0, 1e9)]
med = lambda x: np.median(x) if len(x) else float('nan')
q = lambda x, f: np.percentile(x, f) if len(x) else float('nan')
frac = lambda x: 100. * np.mean(x) if len(x) else float('nan')


def table(title, header, fn):
    print('\n' + title)
    print('  %-9s %-4s %6s ' % ('p [GeV]', 'reg', 'n') + ' '.join('%11s' % h for h in header))
    for reg in ['BTL', 'ETL']:
        for lo, hi in bins:
            s = (region == reg) & (p >= lo) & (p < hi)
            if s.sum() < 20:
                continue
            print('  %-9s %-4s %6d ' % ('%.1f-%s' % (lo, '%.1f' % hi if hi < 1e8 else 'inf'), reg, s.sum()) +
                  ' '.join('%11.4g' % v for v in fn(s)))


print('tracks with an MTD time in both legacy and Alpaka:', len(d), ' matched to truth:', int((ptrue > 0).sum()))

table('1) What the refit does to the momentum (legacy refitted track vs original track, at the PCA)',
      ['med sp/p %', 'med spref/sp', '|dp/p|>10%', '|dp/p|>50%'],
      lambda s: [100 * med(sp[s] / p[s]), med(spref[s] / sp[s]),
                 frac(np.abs(pref[s] / p[s] - 1) > 0.1), frac(np.abs(pref[s] / p[s] - 1) > 0.5)])

t = ptrue > 0
ro = lambda s: np.abs(p[s & t] / ptrue[s & t] - 1)
rr = lambda s: np.abs(pref[s & t] / ptrue[s & t] - 1)
table('2) Is the refitted momentum closer to the truth? |p/p_true - 1| in %, original vs refit (truth-matched)',
      ['q68 orig', 'q68 refit', 'q95 orig', 'q95 refit', '>10% orig', '>10% refit'],
      lambda s: [100 * q(ro(s), 68), 100 * q(rr(s), 68), 100 * q(ro(s), 95), 100 * q(rr(s), 95),
                 frac(ro(s) > 0.1), frac(rr(s) > 0.1)])
table('2b) Pulls (p - p_true)/sigma(p): half width of the central 68%, original vs refit',
      ['pull orig', 'pull refit'],
      lambda s: [0.5 * (q((p - ptrue)[s & t] / sp[s & t], 84) - q((p - ptrue)[s & t] / sp[s & t], 16)),
                 0.5 * (q((pref - ptrue)[s & t] / spref[s & t], 84) - q((pref - ptrue)[s & t] / spref[s & t], 16))])

L = lambda n: col('L_generalTrack' + n)
A = lambda n: col('A_generalTrack' + n)
dpath = lambda s: np.abs(A('PathLength')[s] - L('PathLength')[s])
table('3) Inputs of the time of flight, Alpaka vs legacy (median of the per-track ratio; path differences)',
      ['p A/L', 'path A/L', 'q68 dL cm', '|dL|>1cm %', '|dL|>5cm %', 'sTofPi A/L', 'sTofP A/L'],
      lambda s: [med(A('p')[s] / L('p')[s]), med(A('PathLength')[s] / L('PathLength')[s]), q(dpath(s), 68),
                 frac(dpath(s) > 1), frac(dpath(s) > 5),
                 med(A('SigmaTofPi')[s] / L('SigmaTofPi')[s]), med(A('SigmaTofP')[s] / L('SigmaTofP')[s])])

# sigma(TOF)_legacy should be L/c * sigma(p) * m^2/(p^2 E): which sigma(p) reproduces it?
mp = 0.938272
c_cm_ns = 29.9792458
pl, path = L('p'), L('PathLength')
k = path / c_cm_ns * mp ** 2 / (pl ** 2 * np.sqrt(pl ** 2 + mp ** 2))
table('4) Which sigma(p) does legacy sigma(TOF) correspond to? (median ratio of the implied sigma(p)); legacy p = refit p?',
      ['impl/sp', 'impl/spref', 'pL==pref %'],
      lambda s: [med(L('SigmaTofP')[s] / k[s] / sp[s]), med(L('SigmaTofP')[s] / k[s] / spref[s]),
                 frac(np.abs(pl[s] / pref[s] - 1) < 1e-5)])
