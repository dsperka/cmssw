#!/usr/bin/env python3
"""Summary tables (vertices per event, primary-vertex track composition, vertex time, time per event)
from two harvested MTD validation files, as JSON and as an HTML fragment for render_dqm_gallery.py --summary-html.

  python3 mtd_summary.py --ref DQM_..._legacy_VALID.root --new DQM_..._alpaka_VALID.root \
      --ref-name "Legacy (refit)" --new-name "Alpaka (no refit)" --ms-ref 1102 --ms-new 132 \
      --note "free text shown under the tables" --json summary.json --html summary.html
"""
import argparse
import html
import json

import ROOT

BASE = 'DQMData/Run 1/MTD/Run summary/Vertices/'


def prof_mean(h):
    """Mean of a TProfile over all x bins (per-event average of the profiled quantity)."""
    n = s = 0.
    for b in range(0, h.GetNbinsX() + 2):
        e = h.GetBinEntries(b)
        n += e
        s += h.GetBinContent(b) * e
    return s / n if n else float('nan')


def numbers(fname):
    f = ROOT.TFile.Open(fname)
    g = lambda n: f.Get(BASE + n)
    r = {x: prof_mean(g('PUvs' + x)) for x in ('Real', 'Fake', 'OtherFake', 'Split')}
    r['rec'] = g('RecVerNumber').GetMean()
    r['sel'] = g('RecSelVerNumber').GetMean()
    r['sim'] = g('SimVerNumber').GetMean()
    for n in ('PUTrackRecLVRelMult', 'PUTrackRecLVRelSumPt', 'PUTrackRecLVRelSumWos', 'FakeTrackRecLVRelMult',
              'SecTrackRecLVRelSumWos'):
        h = g(n)
        r[n] = h.GetMean() if h else float('nan')
    t = g('TimeRes')
    r['tRMS'] = t.GetRMS()
    b1, b2 = t.FindBin(-0.03 + 1e-9), t.FindBin(0.03 - 1e-9)
    r['t30'] = t.Integral(b1, b2) / t.Integral(0, t.GetNbinsX() + 1)
    r['tN'] = t.GetEntries()
    ts = g('TimeSignalRes')
    r['tSigRMS'] = ts.GetRMS() if ts else float('nan')
    f.Close()
    return r


def fragment(L, N, a):
    def row(lab, k, fmt='%.3f'):
        return f'<tr><td>{lab}</td><td>{fmt % L[k]}</td><td>{fmt % N[k]}</td></tr>'

    def card(title, head, rows, note=''):
        p = f'<p>{note}</p>' if note else ''
        return (f'<div class="sum-card"><h3>{title}</h3><div class="tbl"><table><tr><th>{head}</th>'
                f'<th>{html.escape(a.ref_name)}</th><th>{html.escape(a.new_name)}</th></tr>{rows}</table></div>{p}</div>')

    verts = (row('Simulated', 'sim', '%.1f') + row('Reco, all', 'rec', '%.1f') + row('Pass ndof', 'sel', '%.1f') +
             row('Real', 'Real', '%.2f') + row('Fake', 'Fake', '%.2f') + row('Split', 'Split', '%.2f') +
             row('Other fake', 'OtherFake', '%.2f') +
             f'<tr><td>Real / (real + fake)</td><td>{L["Real"] / (L["Real"] + L["Fake"]):.3f}</td>'
             f'<td>{N["Real"] / (N["Real"] + N["Fake"]):.3f}</td></tr>'
             f'<tr><td>Real / simulated</td><td>{L["Real"] / L["sim"]:.3f}</td><td>{N["Real"] / N["sim"]:.3f}</td></tr>')
    pv = (row('Pileup tracks, by multiplicity', 'PUTrackRecLVRelMult', '%.4f') +
          row('Pileup tracks, sum pT', 'PUTrackRecLVRelSumPt', '%.4f') +
          row('Pileup tracks, weighted', 'PUTrackRecLVRelSumWos', '%.4f') +
          row('Fake tracks, by multiplicity', 'FakeTrackRecLVRelMult', '%.4f') +
          row('Secondary tracks, weighted', 'SecTrackRecLVRelSumWos', '%.4f'))
    L2 = dict(L, tRMSps=L['tRMS'] * 1e3, tSigps=L['tSigRMS'] * 1e3)
    N2 = dict(N, tRMSps=N['tRMS'] * 1e3, tSigps=N['tSigRMS'] * 1e3)
    tr = lambda lab, k, fmt: f'<tr><td>{lab}</td><td>{fmt % L2[k]}</td><td>{fmt % N2[k]}</td></tr>'
    vt = (tr('RMS of t(rec) - t(sim), all [ps]', 'tRMSps', '%.1f') + tr('RMS, signal vertex [ps]', 'tSigps', '%.1f') +
          tr('Fraction within 30 ps', 't30', '%.3f') + tr('Vertices with a time', 'tN', '%.0f'))
    cards = [card('Vertices per event', 'Per event', verts),
             card('Primary vertex track composition', 'Mean fraction', pv, 'Fractions of the leading-vertex tracks. Lower is better.'),
             card('Vertex time', '', vt)]
    if a.ms_ref is not None and a.ms_new is not None:
        cards.append(card('Time per event', 'ms / event',
                          f'<tr><td>MTD timing module</td><td>{a.ms_ref:.0f}</td><td>{a.ms_new:.0f}</td></tr>',
                          html.escape(a.ms_note) if a.ms_note else ''))
    note = f'<p class="sum-note">{a.note}</p>' if a.note else ''
    return (f'<section class="summary"><h2>Summary{" &middot; " + html.escape(a.events) if a.events else ""}</h2>'
            f'<div class="sum-grid">{"".join(cards)}</div>{note}</section>')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--ref', required=True, help='harvested DQM file, reference (legacy)')
    p.add_argument('--new', required=True, help='harvested DQM file, new (Alpaka)')
    p.add_argument('--ref-name', default='Legacy (refit)')
    p.add_argument('--new-name', default='Alpaka (no refit)')
    p.add_argument('--events', default='', help='e.g. "1000 events"')
    p.add_argument('--ms-ref', type=float)
    p.add_argument('--ms-new', type=float)
    p.add_argument('--ms-note', default='')
    p.add_argument('--note', default='', help='HTML shown under the tables')
    p.add_argument('--json', default='summary.json')
    p.add_argument('--html', default='summary.html')
    a = p.parse_args()
    L, N = numbers(a.ref), numbers(a.new)
    json.dump({'ref': L, 'new': N}, open(a.json, 'w'), indent=1)
    open(a.html, 'w').write(fragment(L, N, a))
    for k in L:
        print('%-24s %12.4f %12.4f' % (k, L[k], N[k]))


if __name__ == '__main__':
    main()
