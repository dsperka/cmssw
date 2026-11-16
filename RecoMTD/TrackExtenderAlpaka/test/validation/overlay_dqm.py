"""Overlay the MTD validation histograms of two harvested DQM files (legacy vs new).
Run inside cmsenv (PyROOT):  python3 overlay_dqm.py legacy.root new.root outdir
Writes one PNG per histogram (1D: overlay + ratio, 2D/profiles: side by side) and summary.csv with
entries, means and a chi2/ndf + KS probability per 1D histogram.
"""
import csv
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

fL, fN, out = ROOT.TFile.Open(sys.argv[1]), ROOT.TFile.Open(sys.argv[2]), sys.argv[3]
os.makedirs(out, exist_ok=True)
base = 'DQMData/Run 1/MTD/Run summary'


def walk(d, path=''):
    for k in d.GetListOfKeys():
        o = k.ReadObj()
        p = path + '/' + k.GetName() if path else k.GetName()
        if o.InheritsFrom('TDirectory'):
            yield from walk(o, p)
        elif o.InheritsFrom('TH1'):
            yield p, o


rows = []
dirL = fL.Get(base)
if not dirL:
    sys.exit('no %s in %s' % (base, sys.argv[1]))
for path, hL in walk(dirL):
    hN = fN.Get(base + '/' + path)
    name = path.replace('/', '__')
    if not hN:
        rows.append([path, 'missing in new', '', '', '', '', '', ''])
        continue
    c = ROOT.TCanvas('c', '', 900, 650)
    if hL.InheritsFrom('TH2') and not hL.InheritsFrom('TProfile'):
        c.Divide(2, 1)
        c.cd(1)
        hL.SetTitle('legacy: ' + hL.GetTitle())
        hL.Draw('colz')
        c.cd(2)
        hN.SetTitle('new: ' + hN.GetTitle())
        hN.Draw('colz')
        rows.append([path, '2D', hL.GetEntries(), hN.GetEntries(), '', '', '', ''])
    else:
        hL.SetLineColor(ROOT.kBlue + 2)
        hN.SetLineColor(ROOT.kOrange + 7)
        hL.SetLineWidth(2)
        hN.SetLineWidth(2)
        pad1 = ROOT.TPad('p1', '', 0, 0.3, 1, 1)
        pad2 = ROOT.TPad('p2', '', 0, 0, 1, 0.3)
        pad1.SetBottomMargin(0.02)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.3)
        pad1.Draw()
        pad2.Draw()
        pad1.cd()
        ymax = max(hL.GetMaximum(), hN.GetMaximum())
        hL.SetMaximum(1.25 * ymax if ymax > 0 else 1)
        hL.Draw('hist' if not hL.InheritsFrom('TProfile') else 'e')
        hN.Draw(('hist' if not hN.InheritsFrom('TProfile') else 'e') + ' same')
        leg = ROOT.TLegend(0.62, 0.78, 0.89, 0.89)
        leg.AddEntry(hL, 'legacy (refit), N=%d' % hL.GetEntries(), 'l')
        leg.AddEntry(hN, 'new (no refit), N=%d' % hN.GetEntries(), 'l')
        leg.SetTextSize(0.035)
        leg.Draw()
        pad2.cd()
        r = hN.Clone('ratio')
        if hN.InheritsFrom('TProfile'):
            r = hN.ProjectionX('ratio')
            r.Divide(hL.ProjectionX('rl'))
        else:
            r.Divide(hL)
        r.SetTitle('')
        r.GetYaxis().SetTitle('new / legacy')
        r.GetYaxis().SetRangeUser(0.8, 1.2)
        r.GetYaxis().SetTitleSize(0.09)
        r.GetYaxis().SetLabelSize(0.08)
        r.GetXaxis().SetLabelSize(0.09)
        r.GetXaxis().SetTitleSize(0.1)
        r.SetLineColor(ROOT.kBlack)
        r.Draw('e')
        chi2 = ks = ''
        if not hL.InheritsFrom('TProfile') and hL.Integral() > 0 and hN.Integral() > 0:
            chi2 = hL.Chi2Test(hN, 'UU CHI2/NDF')
            ks = hL.KolmogorovTest(hN)
        rows.append([path, 'profile' if hL.InheritsFrom('TProfile') else '1D', hL.GetEntries(), hN.GetEntries(),
                     hL.GetMean(), hN.GetMean(), chi2, ks])
    c.SaveAs(os.path.join(out, name + '.png'))
    c.Close()

with open(os.path.join(out, 'summary.csv'), 'w') as f:
    w = csv.writer(f)
    w.writerow(['histogram', 'type', 'entries_legacy', 'entries_new', 'mean_legacy', 'mean_new', 'chi2ndf', 'ks'])
    w.writerows(rows)
print(len(rows), 'histograms')
