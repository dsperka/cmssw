#!/usr/bin/env python3
"""Per-track table for the refit study: legacy (with refit) vs Alpaka (no refit) vs truth.

Per track with an MTD time in both legacy (skim, process RERECO) and Alpaka (chain_cfg.py output, CHAIN):
  original generalTrack p, sigma(p) at the PCA; legacy refitted track (trackExtenderWithMTD, same index) p, sigma(p);
  truth from the closest TrackingParticle (charge, dR < 0.005, |dpt/pt| < 0.3): p, pdg id, production vertex
  time [ns], rho and z, bunch crossing and event id;
  legacy and Alpaka maps: p used for the TOF, path length, sigma(TOF), MTD hit time, t0 and errors.
Writes a flat table (numpy .npz) for analyze_refit.py and analyze_truth.py.
usage: python3 make_table.py skim.root chain.root table.npz
"""
import math
import sys

import numpy as np
import ROOT
from DataFormats.FWLite import Events, Handle

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libFWCoreFWLite')
ROOT.FWLiteEnabler.enable()

MAPS = ['generalTrackp', 'generalTrackPathLength', 'generalTrackSigmaTofPi', 'generalTrackSigmaTofK',
        'generalTrackSigmaTofP', 'generalTracksigmatmtd', 'generalTrackt0', 'generalTracksigmat0',
        'generalTrackTofPi', 'generalTrackTofP', 'generalTracktmtd']


def sigma_p(trk):
    return math.sqrt(max(trk.covariance(0, 0), 0.)) * trk.p() ** 2


def values(vmap):
    # one product (generalTracks) per map: its values in track order
    assert vmap.idSize() == 1
    it = vmap.begin()
    return [it[i] for i in range(it.size())]


def read_maps(ev, proc, handles):
    out = {}
    for name in MAPS:
        ev.getByLabel(('trackExtenderWithMTD', name, proc), handles[name])
        out[name] = handles[name].product()
    return out


def main(skim, chain, outname):
    hT, hR, hTP = Handle('std::vector<reco::io_v1::Track>'), Handle('std::vector<reco::io_v1::Track>'), Handle('std::vector<io_v1::TrackingParticle>')
    hL = {n: Handle('edm::ValueMap<float>') for n in MAPS}
    hA = {n: Handle('edm::ValueMap<float>') for n in MAPS}
    # Alpaka maps, by event number
    alp = {}
    for ev in Events(chain):
        m = read_maps(ev, 'CHAIN', hA)
        alp[ev.eventAuxiliary().event()] = {n: values(m[n]) for n in MAPS}
    rows = []
    for ev in Events(skim):
        e = ev.eventAuxiliary().event()
        if e not in alp:
            continue
        ev.getByLabel(('generalTracks', '', 'RERECO'), hT)
        ev.getByLabel(('trackExtenderWithMTD', '', 'RERECO'), hR)
        ev.getByLabel(('mix', 'MergedTrackTruth', 'HLT'), hTP)
        tracks, refit, tps = hT.product(), hR.product(), hTP.product()
        if refit.size() != tracks.size():
            print('event %d: %d tracks, %d refitted' % (e, tracks.size(), refit.size()))
            continue
        L = {n: values(v) for n, v in read_maps(ev, 'RERECO', hL).items()}
        A = alp[e]
        # truth: charged TPs, binned in eta for a fast nearest search
        tpl = []
        for tp in tps:
            if tp.charge() == 0 or tp.pt() <= 0.3:
                continue
            pv = tp.parentVertex().position()
            tpl.append((tp.eta(), tp.phi(), tp.pt(), tp.p(), tp.charge(),
                        abs(tp.g4Tracks()[0].type()) if tp.g4Tracks().size() else 0,
                        pv.t() * 1e9, math.hypot(pv.x(), pv.y()), pv.z(), tp.eventId().bunchCrossing(), tp.eventId().event()))
        tpa = np.array([t[:5] for t in tpl]) if tpl else np.zeros((0, 5))
        extra = np.array([t[5:] for t in tpl]) if tpl else np.zeros((0, 6))
        for i in range(tracks.size()):
            if L['generalTracksigmatmtd'][i] <= 0 or A['generalTracksigmatmtd'][i] <= 0:
                continue
            t, r = tracks[i], refit[i]
            ptrue = -1
            ex = [-1, -99, -1, -99, -99, -99]
            if len(tpa):
                dphi = np.abs((tpa[:, 1] - t.phi() + math.pi) % (2 * math.pi) - math.pi)
                dr = np.hypot(tpa[:, 0] - t.eta(), dphi)
                ok = (dr < 0.005) & (tpa[:, 4] == t.charge()) & (np.abs(tpa[:, 2] / t.pt() - 1) < 0.3)
                if ok.any():
                    j = np.argmin(np.where(ok, dr, 9))
                    ptrue, ex = tpa[j, 3], list(extra[j])
            rows.append([e, i, t.eta(), t.pt(), t.p(), sigma_p(t), r.p(), sigma_p(r), ptrue] + list(ex) + [t.vz(),
                         t.hitPattern().numberOfValidHits()] +
                        [L[n][i] for n in MAPS] + [A[n][i] for n in MAPS])
    cols = ['event', 'itrk', 'eta', 'pt', 'p', 'sp', 'pref', 'spref', 'ptrue', 'pdg', 'tsim', 'rhosim', 'zsim', 'bx', 'evid', 'vz', 'nhits'] + \
           ['L_' + n for n in MAPS] + ['A_' + n for n in MAPS]
    a = np.array(rows, dtype=float)
    np.savez(outname, cols=np.array(cols), data=a)
    print('tracks:', len(a))


if __name__ == '__main__':
    main(*sys.argv[1:4])
