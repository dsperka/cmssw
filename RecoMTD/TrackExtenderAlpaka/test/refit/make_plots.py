#!/usr/bin/env python3
"""Plots of the refit study (with vs without refit), from the table of make_table.py.

Two questions:
  A. Is the refitted momentum error consistent with the refit itself? (no truth needed)
     If the refit only adds the MTD hit to the original fit (a Kalman update), then
       Var(p_refit - p_orig) = sigma_orig^2 - sigma_refit^2   and   sigma_refit <= sigma_orig,
     so u = (p_refit - p_orig) / sqrt(sigma_orig^2 - sigma_refit^2) is a unit Gaussian.
  B. Against the truth: residuals, pulls, coverage, quoted vs observed resolution; path length and t0.
usage: python3 make_plots.py table.npz outdir
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402
from math import erf, sqrt  # noqa: E402

z = np.load(sys.argv[1])
out = sys.argv[2] if len(sys.argv) > 2 else 'plots'
os.makedirs(out, exist_ok=True)
d = z['data']
c = {n: d[:, i] for i, n in enumerate(z['cols'])}
p, sp, pref, spref, ptrue = c['p'], c['sp'], c['pref'], c['spref'], c['ptrue']
btl = np.abs(c['eta']) < 1.5
REG = [('BTL', btl), ('ETL', ~btl)]
base = p >= 0.7
truth = base & (ptrue > 0)
pion = truth & (c['pdg'] == 211) & (c['bx'] == 0) & (c['rhosim'] >= 0) & (c['rhosim'] < 0.1)
ORIG, REFIT, ALP, GAUS = '#3b6fb6', '#d1495b', '#2a9d8f', '#555555'
PB = np.array([0.7, 1., 1.5, 2., 3., 5., 10.])
plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.3, 'figure.dpi': 110})


def q68(x):
    return 0.5 * (np.percentile(x, 84) - np.percentile(x, 16)) if len(x) > 10 else np.nan


def save(fig, name):
    for ax in fig.axes:
        if ax.get_xscale() == 'log':
            ax.set_xticks([0.7, 1, 1.5, 2, 3, 5, 7])
            ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%g'))
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.tight_layout()
    fig.savefig(os.path.join(out, name + '.png'))
    plt.close(fig)


def gauss(ax, bins, n, width=1.):
    x = 0.5 * (bins[1:] + bins[:-1])
    ax.plot(x, n * (bins[1] - bins[0]) / (sqrt(2 * np.pi) * width) * np.exp(-0.5 * (x / width) ** 2), color=GAUS,
            ls='--', label='unit Gaussian')


def hist2(name, title, xlabel, vals, labels, colors, bins, sel_fn, logy=True, unit_gauss=False):
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for ax, (reg, rm) in zip(axs, REG):
        for v, lab, col in zip(vals, labels, colors):
            m = sel_fn(rm) & np.isfinite(v)
            x = np.clip(v[m], bins[0], bins[-1] - 1e-9)
            ax.hist(x, bins=bins, histtype='step', lw=1.6, color=col,
                    label='%s (q68 %.3g)' % (lab, q68(v[m])))
        if unit_gauss:
            gauss(ax, bins, m.sum())
        ax.set_title('%s: %s' % (reg, title))
        ax.set_xlabel(xlabel)
        if logy:
            ax.set_yscale('log')
            ax.set_ylim(bottom=0.5)
        ax.legend(fontsize=8)
    save(fig, name)


# ---------------- A. consistency of the refit with its own errors ----------------
r = spref / sp
hist2('A1_sigma_ratio', 'refitted / original sigma(p)', r'$\sigma_{refit}(p)/\sigma_{orig}(p)$', [r], ['tracks'],
      [REFIT], np.linspace(0.5, 1.5, 101), lambda rm: base & rm, logy=False)

v = sp ** 2 - spref ** 2
u = np.where(v > 0, (pref - p) / np.sqrt(np.abs(v)), np.nan)
hist2('A2_update_pull', 'refit shift / allowed shift',
      r'$(p_{refit}-p_{orig})/\sqrt{\sigma_{orig}^2-\sigma_{refit}^2}$', [u], ['tracks with sigma shrinking'], [REFIT],
      np.linspace(-8, 8, 81), lambda rm: base & rm & (v > 0), unit_gauss=True)

w = (pref - p) / sp
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for ax, (reg, rm) in zip(axs, REG):
    m = base & rm
    allowed = np.sqrt(np.clip(1 - r[m] ** 2, 0, None))
    bins = np.linspace(0, 5, 101)
    ax.hist(np.clip(np.abs(w[m]), 0, 5 - 1e-9), bins=bins, histtype='step', lw=1.6, color=REFIT,
            label='observed |shift| / sigma_orig')
    # expected: |N(0,1)| * allowed, sampled per track
    rng = np.random.default_rng(1)
    ax.hist(np.clip(np.abs(rng.standard_normal(m.sum())) * allowed, 0, 5 - 1e-9), bins=bins, histtype='step', lw=1.6,
            color=GAUS, ls='--', label='expected from the errors')
    ax.set_yscale('log')
    ax.set_title('%s: how far the refit moves p' % reg)
    ax.set_xlabel(r'$|p_{refit}-p_{orig}|/\sigma_{orig}(p)$')
    ax.legend(fontsize=8)
save(fig, 'A3_shift_vs_expected')

hist2('A4_refit_shift', 'change of p by the refit', r'$(p_{refit}-p_{orig})/p_{orig}$', [(pref - p) / p],
      ['tracks'], [REFIT], np.linspace(-0.5, 0.5, 101), lambda rm: base & rm)

# ---------------- B. against the truth ----------------
ro, rr = p / ptrue - 1, pref / ptrue - 1
hist2('B1_residual', 'momentum residual (prompt pions)', r'$p/p_{true}-1$', [ro, rr], ['original', 'refit'],
      [ORIG, REFIT], np.linspace(-0.3, 0.3, 121), lambda rm: pion & rm)

po, pr = (p - ptrue) / sp, (pref - ptrue) / spref
hist2('B2_pull', 'momentum pull (prompt pions)', r'$(p-p_{true})/\sigma(p)$', [po, pr], ['original', 'refit'],
      [ORIG, REFIT], np.linspace(-10, 10, 101), lambda rm: pion & rm, unit_gauss=True)

fig, axs = plt.subplots(1, 2, figsize=(10, 4))
k = np.linspace(0, 5, 101)
for ax, (reg, rm) in zip(axs, REG):
    m = pion & rm
    for pull, lab, col in ((po, 'original', ORIG), (pr, 'refit', REFIT)):
        a = np.sort(np.abs(pull[m]))
        ax.plot(k, np.searchsorted(a, k) / len(a), color=col, lw=1.8,
                label='%s: %.0f%% within 1 sigma' % (lab, 100 * np.mean(a < 1)))
    ax.plot(k, [erf(x / sqrt(2)) for x in k], color=GAUS, ls='--', label='correct errors: 68% within 1 sigma')
    ax.set_title('%s: coverage of the quoted error (prompt pions)' % reg)
    ax.set_xlabel(r'k: $|p-p_{true}| < k\,\sigma(p)$')
    ax.set_ylabel('fraction of tracks')
    ax.legend(fontsize=8, loc='lower right')
save(fig, 'B3_coverage')

fig, axs = plt.subplots(1, 2, figsize=(10, 4))
xc = np.sqrt(PB[1:] * PB[:-1])
for ax, (reg, rm) in zip(axs, REG):
    for res, sig, lab, col in ((ro, sp, 'original', ORIG), (rr, spref, 'refit', REFIT)):
        obs, quo = [], []
        for lo, hi in zip(PB[:-1], PB[1:]):
            m = pion & rm & (p >= lo) & (p < hi)
            obs.append(100 * q68(res[m]))
            quo.append(100 * np.median(sig[m] / p[m]) if m.sum() else np.nan)
        ax.plot(xc, obs, 'o-', color=col, label='%s: observed (q68 of residual)' % lab)
        ax.plot(xc, quo, 's--', color=col, mfc='none', label='%s: quoted (median sigma(p)/p)' % lab)
    ax.set_xscale('log')
    ax.set_title('%s: momentum resolution, observed vs quoted' % reg)
    ax.set_xlabel('p [GeV]')
    ax.set_ylabel('%')
    ax.legend(fontsize=7)
save(fig, 'B4_resolution_vs_p')

fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for ax, (reg, rm) in zip(axs, REG):
    for pull, lab, col in ((po, 'original', ORIG), (pr, 'refit', REFIT)):
        ax.plot(xc, [q68(pull[pion & rm & (p >= lo) & (p < hi)]) for lo, hi in zip(PB[:-1], PB[1:])], 'o-', color=col,
                label=lab)
    ax.axhline(1, color=GAUS, ls='--', label='correct errors')
    ax.set_xscale('log')
    ax.set_title('%s: pull width (q68) vs p, prompt pions' % reg)
    ax.set_xlabel('p [GeV]')
    ax.legend(fontsize=8)
save(fig, 'B5_pull_width_vs_p')

# ---------------- C. path length and t0: legacy (refit) vs Alpaka (no refit) ----------------
dL = c['A_generalTrackPathLength'] - c['L_generalTrackPathLength']
hist2('C1_path_difference', 'path length, Alpaka - legacy', r'$L_{Alpaka}-L_{legacy}$ [cm]', [dL], ['tracks'],
      [ALP], np.linspace(-40, 40, 161), lambda rm: base & rm)

fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for ax, (reg, rm) in zip(axs, REG):
    for cut, ls in ((1, '-'), (5, '--')):
        ax.plot(xc, [100 * np.mean(np.abs(dL[base & rm & (p >= lo) & (p < hi)]) > cut) for lo, hi in zip(PB[:-1], PB[1:])],
                'o' + ls, color=ALP, label='|dL| > %d cm' % cut)
    ax.set_xscale('log')
    ax.set_title('%s: tracks with different path lengths' % reg)
    ax.set_xlabel('p [GeV]')
    ax.set_ylabel('% of tracks')
    ax.legend(fontsize=8)
save(fig, 'C2_path_difference_vs_p')

same = np.abs(c['L_generalTracktmtd'] - c['A_generalTracktmtd']) < 0.001
tl = 1000 * (c['L_generalTrackt0'] - c['tsim'])
ta = 1000 * (c['A_generalTrackt0'] - c['tsim'])
off = {reg: np.median(tl[pion & rm & same]) for reg, rm in REG}
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for ax, (title, m) in zip(axs, (('BTL, p 0.7-1 GeV, all prompt pions', pion & btl & (p < 1)),
                                ('BTL, same hit, |dL| > 5 cm', pion & btl & same & (np.abs(dL) > 5)))):
    bins = np.linspace(-2000, 2000, 101)
    for t, lab, col in ((tl, 'legacy (refit)', REFIT), (ta, 'Alpaka (no refit)', ALP)):
        x = t[m] - off['BTL']
        ax.hist(np.clip(x, bins[0], bins[-1] - 1e-9), bins=bins, histtype='step', lw=1.6, color=col,
                label='%s: q68 %.0f ps' % (lab, q68(x)))
    ax.set_yscale('log')
    ax.set_title(title)
    ax.set_xlabel('t0 - t_true - %.0f ps (BTL offset) [ps]' % off['BTL'])
    ax.legend(fontsize=8)
save(fig, 'C3_t0_vs_truth_btl')

fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for ax, (reg, rm) in zip(axs, REG):
    for t, lab, col in ((tl, 'legacy (refit)', REFIT), (ta, 'Alpaka (no refit)', ALP)):
        ax.plot(xc, [q68(t[pion & rm & (p >= lo) & (p < hi)]) for lo, hi in zip(PB[:-1], PB[1:])], 'o-', color=col,
                label=lab)
    ax.set_xscale('log')
    ax.set_title('%s: width (q68) of t0 - t_true, prompt pions' % reg)
    ax.set_xlabel('p [GeV]')
    ax.set_ylabel('ps')
    ax.legend(fontsize=8)
save(fig, 'C4_t0_width_vs_p')

fig, ax = plt.subplots(figsize=(6, 4))
bins = np.linspace(-1000, 2000, 151)
for reg, rm in REG:
    x = tl[pion & rm & same]
    ax.hist(np.clip(x, bins[0], bins[-1] - 1e-9), bins=bins, histtype='step', lw=1.6, density=True,
            label='%s: median %+.0f ps' % (reg, np.median(x)))
ax.set_title('t0 - t_true, prompt pions (legacy; Alpaka the same)')
ax.set_xlabel('ps')
ax.legend(fontsize=8)
save(fig, 'C5_btl_offset')

print('plots in', out, ':', ' '.join(sorted(os.listdir(out))))
