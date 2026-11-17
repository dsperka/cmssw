#!/usr/bin/env python3
"""Compare pairs of value maps inside one chain_cfg.py output (device result vs host result).

usage: python3 compare_pairs.py file.root
"""
import sys

import ROOT

from compare_chain import events_of, leaf_values

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

PAIRS = [('floatedmValueMap_trackExtenderWithMTD_mtdQualMVA_CHAIN.', 'floatedmValueMap_mtdTrackQualityMVA_mtdQualMVA_CHAIN.')]
PAIRS += [('floatedmValueMap_tofPID4DnoPIDDevice_%s_CHAIN.' % n, 'floatedmValueMap_tofPID4DnoPID_%s_CHAIN.' % n)
          for n in ['t0', 'sigmat0', 't0safe', 'sigmat0safe', 'probPi', 'probK', 'probP']]


def main():
    f = ROOT.TFile.Open(sys.argv[1])
    t = f.Get('Events')
    t.SetEstimate(20000000)
    ev = events_of(t)
    bad_total = 0
    for a, b in PAIRS:
        if not t.GetBranch(a) or not t.GetBranch(b):
            print('missing %s or %s' % (a, b))
            bad_total += 1
            continue
        va, vb = leaf_values(t, a + 'obj.values_', ev), leaf_values(t, b + 'obj.values_', ev)
        n = bad = 0
        maxdiff = 0.
        examples = []
        for e in va:
            for i, (x, y) in enumerate(zip(va[e], vb[e])):
                n += 1
                if x != y and not (x != x and y != y):
                    bad += 1
                    maxdiff = max(maxdiff, abs(x - y))
                    if len(examples) < 4:
                        examples.append('ev %d trk %d: %.9g vs %.9g' % (e, i, x, y))
            if len(va[e]) != len(vb[e]):
                bad += 1
                examples.append('ev %d: sizes %d vs %d' % (e, len(va[e]), len(vb[e])))
        bad_total += bad
        print('%-60s values %8d differ %6d maxdiff %.3g  %s' % (a.split('_')[1] + ':' + a.split('_')[2], n, bad, maxdiff,
                                                               '; '.join(examples)))
    print('TOTAL differing values: %d' % bad_total)


if __name__ == '__main__':
    main()
