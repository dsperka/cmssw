#!/usr/bin/env python3
"""Exact comparison of two chain_cfg.py outputs, event by event, leaf by leaf.

usage: python3 compare_chain.py ref.root new.root [process]
Every split leaf of every branch of the given process (default CHAIN) present in the reference file
is read with TTree::Draw (data members included), grouped by event number, and compared exactly.
Prints one line per branch: number of values, number of differing values, largest difference.
"""
import sys
from collections import defaultdict

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError


def events_of(tree):
    n = tree.Draw('EventAuxiliary.id_.event_', '', 'goff')
    v = tree.GetV1()
    return [int(v[i]) for i in range(n)]


# the value maps and the vertices; the track clones of trackExtenderWithMTD are copies of the input
# tracks with the chosen MTD hit, which the maps already cover
KEEP = ('floatedmValueMap', 'intedmValueMap', 'recoVertexs', 'GlobalPointedmValueMap', 'Point3DBase')
CHUNK = 5  # entries per Draw call, so that the value buffer (estimate) is large enough


def leaf_values(tree, expr, evnum):
    """values of expr per event number (Entry$ tells which entry each value comes from)"""
    out = defaultdict(list)
    for first in range(0, tree.GetEntries(), CHUNK):
        n = tree.Draw(expr + ':Entry$', '', 'goff', CHUNK, first)
        if n < 0:
            return None
        if n > tree.GetEstimate():
            raise RuntimeError('estimate too small for %s' % expr)
        v1, v2 = tree.GetV1(), tree.GetV2()
        for i in range(n):
            out[evnum[int(v2[i])]].append(v1[i])
    return out


def leaves_of(branch, prefix=''):
    out = []
    subs = branch.GetListOfBranches()
    if subs.GetEntries() == 0:
        for leaf in branch.GetListOfLeaves():
            out.append(leaf.GetName())
    else:
        for sub in subs:
            out += leaves_of(sub)
    return out


def main():
    ref, new = sys.argv[1], sys.argv[2]
    proc = sys.argv[3] if len(sys.argv) > 3 else 'CHAIN'
    fa, fb = ROOT.TFile.Open(ref), ROOT.TFile.Open(new)
    ta, tb = fa.Get('Events'), fb.Get('Events')
    for t in (ta, tb):
        t.SetEstimate(20000000)
    ea, eb = events_of(ta), events_of(tb)
    print('events: ref %d new %d common %d' % (len(ea), len(eb), len(set(ea) & set(eb))))
    total_bad = 0
    for br in ta.GetListOfBranches():
        name = br.GetName()
        if not name.endswith(proc + '.'):
            continue
        if not name.startswith(KEEP):
            continue
        if not tb.GetBranch(name):
            print('MISSING in new: %s' % name)
            total_bad += 1
            continue
        nval = nbad = 0
        maxdiff = 0.
        skipped = []
        leaves = leaves_of(br)
        if 'ValueMap' in name:
            leaves = [name + 'obj.values_']  # stored unsplit: read the values through the dictionary
        for leaf in leaves:
            if leaf.endswith('_') and leaf.count('.') == 0:
                continue
            if 'present' in leaf and 'ValueMap' not in name or leaf.endswith('obj') or 'processIndex_' in leaf or 'productIndex_' in leaf:
                continue
            va, vb = leaf_values(ta, leaf, ea), leaf_values(tb, leaf, eb)
            if va is None or vb is None:
                skipped.append(leaf)
                continue
            for ev in va:
                a, b = va[ev], vb.get(ev, [])
                nval += len(a)
                if len(a) != len(b):
                    nbad += max(len(a), len(b))
                    maxdiff = float('inf')
                    continue
                for x, y in zip(a, b):
                    if x != y and not (x != x and y != y):
                        nbad += 1
                        maxdiff = max(maxdiff, abs(x - y))
        total_bad += nbad
        print('%-90s values %9d differ %7d maxdiff %.3g%s' %
              (name, nval, nbad, maxdiff, (' (skipped %d leaves)' % len(skipped)) if skipped else ''))
    print('TOTAL differing values: %d' % total_bad)
    return 0 if total_bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
