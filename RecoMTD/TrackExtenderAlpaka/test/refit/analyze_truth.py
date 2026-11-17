#!/usr/bin/env python3
"""Track t0 (pion hypothesis) of legacy and Alpaka against the true production time.

Sample: truth-matched charged pions from in-time interactions (bx 0) produced near the beam line
(rho < 0.1 cm), with an MTD time in both legacy and Alpaka.
Split into tracks where legacy and Alpaka use the same MTD hit (|tmtd diff| < 1 ps) or not.
usage: python3 analyze_truth.py table.npz      (table from make_table.py)
"""
import sys

import numpy as np

z = np.load(sys.argv[1])
cols = list(z['cols'])
d = z['data']
c = {n: d[:, i] for i, n in enumerate(cols)}

sel = (c['pdg'] == 211) & (c['bx'] == 0) & (c['rhosim'] < 0.1) & (c['rhosim'] >= 0)
print('tracks with MTD time in both: %d, truth-matched: %d, prompt in-time pions: %d'
      % (len(d), (c['ptrue'] > 0).sum(), sel.sum()))
btl = np.abs(c['eta']) < 1.5
same = np.abs(c['L_generalTracktmtd'] - c['A_generalTracktmtd']) < 0.001
dL = c['A_generalTrackPathLength'] - c['L_generalTrackPathLength']
rL = c['L_generalTrackt0'] - c['tsim']
rA = c['A_generalTrackt0'] - c['tsim']
pL = rL / c['L_generalTracksigmat0']
pA = rA / c['A_generalTracksigmat0']
print('same MTD hit: %.1f%% of the selected pions' % (100. * same[sel].mean()))


def q68(x):
    return 0.5 * (np.percentile(x, 84) - np.percentile(x, 16)) if len(x) else np.nan


def line(name, m):
    n = m.sum()
    if n < 20:
        return
    out = []
    for r, p in ((rL, pL), (rA, pA)):
        x = r[m] * 1000.
        out.append('%6.1f %6.1f %6.1f %5.2f %5.1f' % (np.median(x), q68(x), np.sqrt(np.mean(x ** 2)),
                                                       q68(p[m]), 100. * np.mean(np.abs(p[m]) > 3)))
    print('  %-14s %7d | %s | %s' % (name, n, out[0], out[1]))


hdr = '  %-14s %7s | %s | %s' % ('', 'tracks', 'med   q68   rms  pull  >3s%', 'med   q68   rms  pull  >3s%')
PB = [(0.7, 1.), (1., 1.5), (1.5, 2.), (2., 3.), (3., 5.), (5., 1e9)]
for reg, rm in (('BTL', btl), ('ETL', ~btl)):
    for label, hm in (('all', np.ones(len(d), bool)), ('same MTD hit', same), ('different hit', ~same)):
        print('\n%s, %s: t0 - t_sim in ps (legacy | Alpaka); pull = q68 of residual/sigma(t0)' % (reg, label))
        print(hdr)
        for lo, hi in PB:
            line('p %.1f-%s' % (lo, '%.1f' % hi if hi < 100 else 'inf'), sel & rm & hm & (c['p'] >= lo) & (c['p'] < hi))
        line('all p', sel & rm & hm & (c['p'] >= 0.7))

print('\nTime offset per region (prompt in-time pions, same MTD hit, p > 0.7 GeV): medians in ps')
off = {}
for reg, rm in (('BTL', btl), ('ETL', ~btl)):
    m = sel & rm & same & (c['p'] >= 0.7)
    off[reg] = np.median(rL[m])
    print('  %s: t0 - t_sim %.0f; t_mtd - t_sim %.0f, TOF(pi) %.0f; sigma(t_mtd) %.0f, sigma(t0) legacy %.0f Alpaka %.0f'
          % (reg, 1000 * off[reg], 1000 * np.median((c['L_generalTracktmtd'] - c['tsim'])[m]),
             1000 * np.median(c['L_generalTrackTofPi'][m]), 1000 * np.median(c['L_generalTracksigmatmtd'][m]),
             1000 * np.median(c['L_generalTracksigmat0'][m]), 1000 * np.median(c['A_generalTracksigmat0'][m])))

print('\nSame MTD hit, path lengths differ: whose t0 is closer to the truth? (region offset above removed)')
for reg, rm in (('BTL', btl), ('ETL', ~btl)):
    for lo, hi in ((1., 5.), (5., 1e9)):
        m = sel & rm & same & (c['p'] >= 0.7) & (np.abs(dL) >= lo) & (np.abs(dL) < hi)
        if m.sum() < 5:
            continue
        aL, aA = np.abs(rL[m] - off[reg]), np.abs(rA[m] - off[reg])
        print('  %s |dL| %g-%s cm: %5d tracks (median p %.2f GeV, median dL %+.1f cm); Alpaka closer %.0f%%; '
              'median |res| legacy %.0f ps, Alpaka %.0f ps'
              % (reg, lo, '%g' % hi if hi < 100 else 'inf', m.sum(), np.median(c['p'][m]), np.median(dL[m]),
                 100. * np.mean(aA < aL), 1000. * np.median(aL), 1000. * np.median(aA)))
