#!/usr/bin/env python3
"""
Self-contained PV DQM gallery generator.

Compares harvested DQM files (from `harvester.py`, via `cmsRun testVertexing.py -d`)
across however many of Old / InBlocks / Alpaka you ran, and writes a single
self-contained HTML page with one overlay plot per histogram common to all
of them.

Usage:
  python3 render_dqm_gallery.py \
      --old Old_harvested.root \
      --inblocks InBlocks_harvested.root \
      --alpaka Alpaka_harvested.root \
      --out pv_dqm_gallery.html

Any subset of --old/--inblocks/--alpaka may be omitted. Collection names
default to "offlinePrimaryVertices" for Old/InBlocks and "vertexAoS" for
Alpaka (override with --old-collection/--inblocks-collection/--alpaka-collection
if your config uses different labels).
"""
import argparse
import base64
import html
import json
import os
import re
import sys
import tempfile

import ROOT

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import dqm_labels  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

COLORS = {
    "Old": (ROOT.kBlack, ROOT.kFullCircle),
    "InBlocks": (ROOT.kRed + 1, ROOT.kFullSquare),
    "Alpaka": (ROOT.kAzure + 2, ROOT.kFullTriangleUp),
}
SERIES_ORDER = ["Old", "InBlocks", "Alpaka"]


def categorize(name):
    if name.startswith(("effic_", "fakerate_", "duplicate_", "merged_", "gen_duplicate_")):
        return "Efficiency, fake, duplicate & merge rates"
    if name.startswith("MisTagRate"):
        return "Mis-tag rate"
    if "_Pull" in name or "_Resol" in name:
        return "Pulls & resolutions"
    if name.startswith(("TruePVLocationIndex", "KindOfSignalPV", "PV_purity_vs_", "globalEfficiencies")):
        return "PV identification"
    if name.startswith(("GenAllAssoc2Reco", "GenAllAssoc2RecoMatched", "GenAllAssoc2RecoMultiMatched", "GenPVAssoc2RecoPV")):
        return "Gen-vertex distributions"
    if name.startswith(("RecoAllAssoc2Gen", "RecoAllAssoc2MultiMatchedGen", "RecoAssoc2GenPV", "RecoPVAssoc2GenPV")):
        return "Reco-vertex distributions"
    return "Other"


CATEGORY_ORDER = [
    "Efficiency, fake, duplicate & merge rates",
    "Mis-tag rate",
    "Pulls & resolutions",
    "Gen-vertex distributions",
    "Reco-vertex distributions",
    "PV identification",
    "Other",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--old", help="Old mode harvested DQM root file")
    p.add_argument("--inblocks", help="InBlocks mode harvested DQM root file")
    p.add_argument("--alpaka", help="Alpaka mode harvested DQM root file")
    p.add_argument("--old-collection", default="offlinePrimaryVertices")
    p.add_argument("--inblocks-collection", default="offlinePrimaryVertices")
    p.add_argument("--alpaka-collection", default="vertexAoS")
    p.add_argument("--out", default="pv_dqm_gallery.html", help="Output HTML file")
    p.add_argument("--title", default="PrimaryVertexProducer DQM comparison")
    p.add_argument("--subtitle", default="", help="Extra text for the page subtitle")
    p.add_argument("--workdir", default=None, help="Directory to write intermediate PNGs (default: temp dir)")
    return p.parse_args()


def build_series(args):
    series = {}
    if args.old:
        series["Old"] = (args.old, args.old_collection)
    if args.inblocks:
        series["InBlocks"] = (args.inblocks, args.inblocks_collection)
    if args.alpaka:
        series["Alpaka"] = (args.alpaka, args.alpaka_collection)
    if len(series) < 2:
        sys.exit("Need at least two of --old/--inblocks/--alpaka to compare.")
    return series


def render(series, outdir, label_order):
    dirs = {}
    tfiles = {}
    for label, (fn, coll) in series.items():
        tf = ROOT.TFile.Open(fn)
        if not tf or tf.IsZombie():
            sys.exit(f"Could not open {fn}")
        tfiles[label] = tf
        d = tf.Get(f"DQMData/Run 1/Vertexing/Run summary/PrimaryVertexV/{coll}")
        if not d:
            sys.exit(f"Could not find PrimaryVertexV/{coll} in {fn}")
        dirs[label] = d

    namesets = {}
    classmap = {}
    for label, d in dirs.items():
        names = set()
        for k in d.GetListOfKeys():
            cn = k.GetClassName()
            if cn in ("TH1F", "TProfile"):
                names.add(k.GetName())
                classmap[k.GetName()] = cn
        namesets[label] = names

    common = sorted(set.intersection(*namesets.values()))
    print(f"Total common 1D histograms: {len(common)}", file=sys.stderr)

    manifest = []
    safe_re = re.compile(r"[^A-Za-z0-9_.-]")

    for i, name in enumerate(common):
        hists = {}
        for label, d in dirs.items():
            h = d.Get(name)
            h.SetDirectory(0)
            hists[label] = h

        first_label = label_order[0]
        nbins = hists[first_label].GetNbinsX()
        first_populated, last_populated = None, None
        for b in range(1, nbins + 1):
            if any(hists[l].GetBinContent(b) != 0 for l in hists):
                if first_populated is None:
                    first_populated = b
                last_populated = b
        if first_populated is None:
            first_populated, last_populated = 1, nbins

        pad = max(1, int(round((last_populated - first_populated) * 0.08)))
        lo_bin = max(1, first_populated - pad)
        hi_bin = min(nbins, last_populated + pad)

        overflow_notes = []
        for label in label_order:
            h = hists[label]
            total = h.GetEntries()
            if total <= 0:
                continue
            outside = h.GetBinContent(0) + h.GetBinContent(nbins + 1)
            frac = outside / total
            if frac > 0.02:
                overflow_notes.append(f"{label} {frac*100:.0f}% outside axis range")

        canvas = ROOT.TCanvas(f"c_{i}", "", 600, 460)
        canvas.SetLeftMargin(0.16)
        canvas.SetBottomMargin(0.16)
        canvas.SetTopMargin(0.15 if overflow_notes else 0.06)
        canvas.SetRightMargin(0.04)
        canvas.SetGrid(1, 1)

        ymax, ymin = 0.0, 0.0
        for label in hists:
            h = hists[label]
            for b in range(lo_bin, hi_bin + 1):
                v = h.GetBinContent(b)
                e = h.GetBinError(b)
                if v + e > ymax:
                    ymax = v + e
                if v - e < ymin:
                    ymin = v - e
        headroom = (ymax - ymin) * 0.15 if ymax > ymin else 1.0

        _title, xlabel, ylabel, _definition = dqm_labels.derive(name)
        if not xlabel:
            xlabel = name

        if overflow_notes:
            leg = ROOT.TLegend(0.55, 0.70, 0.95, 0.88)
        else:
            leg = ROOT.TLegend(0.55, 0.78, 0.95, 0.94)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.042)

        first = True
        for label in label_order:
            h = hists[label]
            color, marker = COLORS[label]
            h.SetLineColor(color)
            h.SetMarkerColor(color)
            h.SetMarkerStyle(marker)
            h.SetMarkerSize(0.75)
            h.SetLineWidth(2)
            h.GetXaxis().SetRangeUser(h.GetXaxis().GetBinLowEdge(lo_bin), h.GetXaxis().GetBinUpEdge(hi_bin))
            h.GetXaxis().SetLabelSize(0.042)
            h.GetYaxis().SetLabelSize(0.042)
            h.GetXaxis().SetTitle(xlabel)
            h.GetYaxis().SetTitle(ylabel)
            h.GetXaxis().SetTitleSize(0.046)
            h.GetYaxis().SetTitleSize(0.046)
            h.GetXaxis().SetTitleOffset(1.35)
            h.GetYaxis().SetTitleOffset(1.55)
            h.SetMaximum(ymax + headroom)
            h.SetMinimum(ymin - headroom if ymin < 0 else 0)
            opt = "E1" if first else "E1 SAME"
            h.Draw(opt)
            leg.AddEntry(h, label, "lp")
            first = False

        leg.Draw()

        if overflow_notes:
            note = ROOT.TLatex()
            note.SetNDC(True)
            note.SetTextFont(42)
            note.SetTextSize(0.032)
            note.SetTextColor(ROOT.kGray + 3)
            for j, line in enumerate(overflow_notes):
                note.DrawLatex(0.17, 0.965 - j * 0.042, line)

        safe_name = safe_re.sub("_", name)
        outpath = os.path.join(outdir, safe_name + ".png")
        canvas.Print(outpath)
        canvas.Close()

        manifest.append({
            "name": name,
            "file": safe_name + ".png",
            "category": categorize(name),
            "class": classmap[name],
        })

        if (i + 1) % 50 == 0:
            print(f"...{i+1}/{len(common)}", file=sys.stderr)

    for tf in tfiles.values():
        tf.Close()

    print(f"DONE {len(manifest)}", file=sys.stderr)
    return manifest


def build_html(manifest, outdir, series, title, subtitle):
    manifest = sorted(manifest, key=lambda x: (x["category"], x["name"]))
    by_cat = {}
    for item in manifest:
        by_cat.setdefault(item["category"], []).append(item)
    cat_counts = {c: len(by_cat.get(c, [])) for c in CATEGORY_ORDER if by_cat.get(c)}

    cards_html = []
    for cat in CATEGORY_ORDER:
        items = by_cat.get(cat)
        if not items:
            continue
        cat_id = cat.lower().replace(" ", "-").replace(",", "").replace("&", "and")
        cards_html.append(f'<section class="cat-section" id="cat-{cat_id}" data-cat="{html.escape(cat)}">')
        cards_html.append(f'<h2 class="cat-heading">{html.escape(cat)} <span class="cat-count">{len(items)}</span></h2>')
        cards_html.append('<div class="grid">')
        for item in items:
            with open(os.path.join(outdir, item["file"]), "rb") as imgf:
                b64 = base64.b64encode(imgf.read()).decode("ascii")
            raw_name = item["name"]
            t, _xlab, _ylab, definition = dqm_labels.derive(raw_name)
            name = html.escape(raw_name)
            title_esc = html.escape(t)
            def_html = f'<div class="card-def">{html.escape(definition)}</div>' if definition else ""
            cards_html.append(
                f'<figure class="card" data-name="{name.lower()} {title_esc.lower()}" data-cat="{html.escape(cat)}">'
                f'<img src="data:image/png;base64,{b64}" alt="{title_esc}" loading="lazy">'
                f'<figcaption>'
                f'<div class="card-title">{title_esc}</div>'
                f'{def_html}'
                f'<div class="card-rawname">{name}</div>'
                f'</figcaption>'
                f'</figure>'
            )
        cards_html.append('</div></section>')

    chips_html = ['<button class="chip active" data-filter="all">All <span>{}</span></button>'.format(len(manifest))]
    for cat in CATEGORY_ORDER:
        if cat_counts.get(cat):
            cat_id = cat.lower().replace(" ", "-").replace(",", "").replace("&", "and")
            chips_html.append(f'<button class="chip" data-filter="{cat_id}">{html.escape(cat)} <span>{cat_counts[cat]}</span></button>')

    legend_items = []
    for label in SERIES_ORDER:
        if label not in series:
            continue
        color_hex = {"Old": "#000", "InBlocks": "#e3281f", "Alpaka": "#2878c8"}[label]
        legend_items.append(f'<span><i class="dot" style="background:{color_hex}"></i> {html.escape(label)}</span>')

    CARDS = "\n".join(cards_html)
    CHIPS = "\n".join(chips_html)
    LEGEND = "\n      ".join(legend_items)
    TOTAL = len(manifest)
    SUBTITLE_EXTRA = f" &middot; {html.escape(subtitle)}" if subtitle else ""

    page = f"""<title>{html.escape(title)}</title>
<style>
  :root {{
    --bg: #f3f1ec; --surface: #ffffff; --ink: #1c1d20; --ink-dim: #5c5e63;
    --rule: #ddd9d0; --accent: #2d6e6a; --accent-ink: #1d4b48; --chip-bg: #eae7df;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{ --bg: #15171a; --surface: #1c1f23; --ink: #e9e7e1; --ink-dim: #9b9d9f;
      --rule: #32353a; --accent: #5bb3ac; --accent-ink: #8ad4cd; --chip-bg: #24272b; }}
  }}
  :root[data-theme="dark"] {{ --bg: #15171a; --surface: #1c1f23; --ink: #e9e7e1; --ink-dim: #9b9d9f;
    --rule: #32353a; --accent: #5bb3ac; --accent-ink: #8ad4cd; --chip-bg: #24272b; }}
  :root[data-theme="light"] {{ --bg: #f3f1ec; --surface: #ffffff; --ink: #1c1d20; --ink-dim: #5c5e63;
    --rule: #ddd9d0; --accent: #2d6e6a; --accent-ink: #1d4b48; --chip-bg: #eae7df; }}

  * {{ box-sizing: border-box; }}
  html, body {{ margin: 0; padding: 0; }}
  body {{ background: var(--bg); color: var(--ink);
    font-family: ui-sans-serif, -apple-system, "Helvetica Neue", Arial, sans-serif; -webkit-font-smoothing: antialiased; }}
  .wrap {{ max-width: 1400px; margin: 0 auto; padding: 0 24px 64px; }}
  header.top {{ position: sticky; top: 0; z-index: 10; background: var(--bg); border-bottom: 1px solid var(--rule); padding: 20px 24px 14px; }}
  .top-inner {{ max-width: 1400px; margin: 0 auto; }}
  .eyebrow {{ font-family: ui-monospace, "SF Mono", "Cascadia Code", monospace; font-size: 0.7rem; letter-spacing: 0.08em;
    text-transform: uppercase; color: var(--accent-ink); margin: 0 0 6px; }}
  h1 {{ font-size: 1.5rem; font-weight: 650; margin: 0 0 4px; letter-spacing: -0.01em; }}
  .sub {{ color: var(--ink-dim); font-size: 0.88rem; margin: 0 0 14px; }}
  .legend-key {{ display: flex; gap: 18px; font-size: 0.8rem; margin-bottom: 14px; flex-wrap: wrap; }}
  .legend-key span {{ display: flex; align-items: center; gap: 6px; }}
  .dot {{ width: 10px; height: 10px; border-radius: 50%; display: inline-block; }}
  .controls {{ display: flex; gap: 10px; align-items: center; flex-wrap: wrap; }}
  #search {{ flex: 1; min-width: 220px; max-width: 360px; padding: 8px 12px; border-radius: 6px; border: 1px solid var(--rule);
    background: var(--surface); color: var(--ink); font-size: 0.88rem; }}
  #search:focus {{ outline: 2px solid var(--accent); outline-offset: 1px; }}
  .chips {{ display: flex; gap: 8px; flex-wrap: wrap; }}
  .chip {{ font-family: ui-monospace, monospace; font-size: 0.72rem; background: var(--chip-bg); border: 1px solid var(--rule);
    color: var(--ink); border-radius: 999px; padding: 6px 12px; cursor: pointer; }}
  .chip span {{ color: var(--ink-dim); margin-left: 4px; }}
  .chip.active {{ background: var(--accent); color: #fff; border-color: var(--accent); }}
  .chip.active span {{ color: rgba(255,255,255,0.8); }}
  .cat-heading {{ font-size: 1.1rem; font-weight: 650; margin: 36px 0 14px; display: flex; align-items: baseline; gap: 10px; }}
  .cat-count {{ font-family: ui-monospace, monospace; font-size: 0.75rem; color: var(--ink-dim); font-weight: 400; }}
  .grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); gap: 14px; }}
  .card {{ margin: 0; background: var(--surface); border: 1px solid var(--rule); border-radius: 8px; padding: 8px; overflow: hidden; }}
  .card img {{ width: 100%; display: block; border-radius: 4px; }}
  .card figcaption {{ margin-top: 8px; }}
  .card-title {{ font-size: 0.86rem; font-weight: 650; color: var(--ink); text-align: center; line-height: 1.3; }}
  .card-def {{ font-size: 0.74rem; color: var(--ink-dim); text-align: center; margin-top: 3px; line-height: 1.35; }}
  .card-rawname {{ font-family: ui-monospace, monospace; font-size: 0.64rem; color: var(--ink-dim); opacity: 0.6; margin-top: 5px;
    word-break: break-word; text-align: center; }}
  .card.hidden {{ display: none; }}
  .cat-section.hidden {{ display: none; }}
  .glossary {{ background: var(--surface); border: 1px solid var(--rule); border-radius: 8px; padding: 14px 18px; margin: 14px 0 4px; font-size: 0.82rem; }}
  .glossary summary {{ cursor: pointer; font-weight: 650; }}
  .glossary dl {{ margin: 12px 0 0; display: grid; grid-template-columns: max-content 1fr; gap: 6px 16px; }}
  .glossary dt {{ font-weight: 650; white-space: nowrap; }}
  .glossary dd {{ margin: 0; color: var(--ink-dim); }}
  footer {{ margin-top: 48px; padding-top: 16px; border-top: 1px solid var(--rule); font-family: ui-monospace, monospace;
    font-size: 0.72rem; color: var(--ink-dim); }}
</style>

<header class="top">
  <div class="top-inner">
    <p class="eyebrow">RecoVertex / PrimaryVertexProducer &middot; DQM comparison</p>
    <h1>{html.escape(title)}</h1>
    <p class="sub">{TOTAL} histograms{SUBTITLE_EXTRA} &middot; axes cropped to the populated data range, with an overflow note where a histogram's fixed binning cuts off part of the distribution</p>
    <div class="legend-key">
      {LEGEND}
    </div>
    <div class="controls">
      <input id="search" type="text" placeholder="Filter by histogram name&hellip;">
      <div class="chips">
        {CHIPS}
      </div>
    </div>
  </div>
</header>

<div class="wrap">
<details class="glossary">
  <summary>What do these metrics mean?</summary>
  <dl>
    <dt>Efficiency</dt><dd>Fraction of truth vertices with a matching reconstructed vertex.</dd>
    <dt>Fake rate</dt><dd>Fraction of reconstructed vertices with no matching truth vertex.</dd>
    <dt>Duplicate rate</dt><dd>Fraction of reconstructed vertices whose truth vertex is also claimed by another reconstructed vertex.</dd>
    <dt>Gen-level duplicate rate</dt><dd>Same split, counted from the truth side.</dd>
    <dt>Merge rate</dt><dd>Fraction of reconstructed vertices that match more than one truth vertex.</dd>
    <dt>Mis-tag rate</dt><dd>Fraction of events where the true signal vertex was not ranked as the primary reconstructed vertex.</dd>
    <dt>Pull</dt><dd>(reconstructed &minus; true) divided by the estimated uncertainty; should look like a unit Gaussian if uncertainties are well calibrated.</dd>
    <dt>Resolution</dt><dd>Reconstructed &minus; true, in the quantity's own units; measures accuracy directly.</dd>
    <dt>"Vs." plots</dt><dd>The rate or resolution metric, binned by the variable named after "vs."</dd>
  </dl>
</details>
{CARDS}
<footer>
  DataFormats/VertexSoA &middot; RecoVertex/PrimaryVertexProducer &middot; generated by render_dqm_gallery.py
</footer>
</div>

<script>
(function() {{
  var search = document.getElementById('search');
  var chips = document.querySelectorAll('.chip');
  var sections = document.querySelectorAll('.cat-section');
  var activeFilter = 'all';
  function apply() {{
    var q = search.value.trim().toLowerCase();
    sections.forEach(function(sec) {{
      var secId = sec.id;
      var matchCat = activeFilter === 'all' || secId === 'cat-' + activeFilter;
      var anyVisible = false;
      sec.querySelectorAll('.card').forEach(function(card) {{
        var matchName = !q || card.getAttribute('data-name').indexOf(q) !== -1;
        var show = matchCat && matchName;
        card.classList.toggle('hidden', !show);
        if (show) anyVisible = true;
      }});
      sec.classList.toggle('hidden', !anyVisible);
    }});
  }}
  search.addEventListener('input', apply);
  chips.forEach(function(chip) {{
    chip.addEventListener('click', function() {{
      chips.forEach(function(c) {{ c.classList.remove('active'); }});
      chip.classList.add('active');
      activeFilter = chip.getAttribute('data-filter');
      apply();
    }});
  }});
}})();
</script>
"""
    return page


def main():
    args = parse_args()
    series = build_series(args)
    label_order = [l for l in SERIES_ORDER if l in series]

    if args.workdir:
        os.makedirs(args.workdir, exist_ok=True)
        manifest = render(series, args.workdir, label_order)
        page = build_html(manifest, args.workdir, series, args.title, args.subtitle)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = render(series, tmp, label_order)
            page = build_html(manifest, tmp, series, args.title, args.subtitle)

    with open(args.out, "w") as f:
        f.write(page)
    print(f"wrote {args.out} ({os.path.getsize(args.out)} bytes)", file=sys.stderr)


if __name__ == "__main__":
    main()
