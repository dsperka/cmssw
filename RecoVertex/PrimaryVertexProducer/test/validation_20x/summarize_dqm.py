"""Per-event real/fake/split vertex averages from a harvested MTD vertex validation file.
   usage: python3 summarize_dqm.py <harvested.root>"""
import ROOT, sys
f=ROOT.TFile(sys.argv[1] if len(sys.argv) > 1 else "DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root")
base="DQMData/Run 1/MTD/Run summary/Vertices"
d=f.Get(base); print([k.GetName() for k in d.GetListOfKeys()])
def prof_mean(h):
    n=0.;s=0.
    for b in range(0,h.GetNbinsX()+2):
        e=h.GetBinEntries(b); n+=e; s+=h.GetBinContent(b)*e
    return s/n if n else float("nan")
cols=["3D","4D","GNN","InBlocks","Alpaka"]
res={}
for c in cols:
    g=d.Get(c)
    if not g: print("MISSING",c); continue
    r={k:prof_mean(g.Get("PUvs"+k)) for k in ("Real","Fake","OtherFake","Split")}
    r["rec"]=g.Get("RecVerNumber").GetMean(); r["sel"]=g.Get("RecSelVerNumber").GetMean(); r["sim"]=g.Get("SimVerNumber").GetMean()
    res[c]=r
print("%-11s"%"", "  ".join("%9s"%c for c in res))
for k,l in (("sim","sim"),("rec","reco all"),("sel","pass ndof"),("Real","real"),("Fake","fake"),("Split","split"),("OtherFake","other fake")):
    print("%-11s"%l,"  ".join("%9.2f"%res[c][k] for c in res))
print("%-11s"%"real/(r+f)","  ".join("%9.3f"%(res[c]["Real"]/(res[c]["Real"]+res[c]["Fake"])) for c in res))
print("%-11s"%"split frac","  ".join("%9.3f"%(res[c]["Split"]/(res[c]["Real"]+res[c]["Fake"])) for c in res))
print("%-11s"%"real/sim","  ".join("%9.3f"%(res[c]["Real"]/res[c]["sim"]) for c in res))
