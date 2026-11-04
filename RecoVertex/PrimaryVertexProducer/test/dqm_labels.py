import re

VAR_LABELS = {
    "NumTracks": "number of tracks",
    "NumVertices": "number of vertices",
    "ClosestDistanceZ": "distance to closest vertex in z",
    "ClosestVertexInZ": "distance to closest vertex in z",
    "PairDistanceZ": "distance between vertex pair in z",
    "X": "vertex x",
    "Y": "vertex y",
    "Z": "vertex z",
    "R": "vertex radial distance",
    "Pt2": "sum p_{T}^{2} of tracks",
    "Pt": "track p_{T}",
    "PU": "number of pileup interactions",
    "Purity": "vertex purity",
    "Chi2": "chi-square",
    "Chi2Ndof": "chi-square / ndof",
    "Ndof": "degrees of freedom",
    "PullX": "x pull",
    "PullY": "y pull",
    "PullZ": "z pull",
    "ResolX": "x resolution",
    "ResolY": "y resolution",
    "ResolZ": "z resolution",
    "ResolPt2": "sum p_{T}^{2} resolution",
    "SharedTrackFractionAssociationScore": "shared track fraction association score",
    "Missing": "missing fraction",
    "Index": "vertex index",
}

AXIS_LABELS = {k: v[0].upper() + v[1:] for k, v in VAR_LABELS.items()}
AXIS_LABELS["ClosestDistanceZ"] = "Distance to closest vertex in z [cm]"
AXIS_LABELS["ClosestVertexInZ"] = "Distance to closest vertex in z [cm]"
AXIS_LABELS["PairDistanceZ"] = "Distance between vertex pair in z [cm]"
AXIS_LABELS["X"] = "Vertex x [cm]"
AXIS_LABELS["Y"] = "Vertex y [cm]"
AXIS_LABELS["Z"] = "Vertex z [cm]"
AXIS_LABELS["R"] = "Vertex radial distance [cm]"
AXIS_LABELS["Pt2"] = "Sum p_{T}^{2} of tracks [GeV^{2}]"
AXIS_LABELS["Pt"] = "Track p_{T} [GeV]"

RATE_DEFS = {
    "effic_vs_": ("Efficiency", "Fraction of truth vertices with a matching reconstructed vertex."),
    "fakerate_vs_": ("Fake rate", "Fraction of reconstructed vertices with no matching truth vertex."),
    "gen_duplicate_vs_": ("Gen-level duplicate rate", "Fraction of truth vertices matched by more than one reconstructed vertex."),
    "duplicate_vs_": ("Duplicate rate", "Fraction of reconstructed vertices whose truth vertex is also claimed by another reconstructed vertex (the truth vertex was split into two reco vertices)."),
    "merged_vs_": ("Merge rate", "Fraction of reconstructed vertices that match more than one truth vertex (several truth vertices merged into one reco vertex)."),
}

PREFIX_LABELS = [
    ("GenAllAssoc2RecoMultiMatched", "Truth vertices matched by more than one reconstructed vertex"),
    ("GenAllAssoc2RecoMatched", "Truth vertices matched to a reconstructed vertex"),
    ("GenAllAssoc2Reco", "All truth vertices"),
    ("GenPVAssoc2RecoPVMatched", "True signal vertex, matched to the reconstructed PV"),
    ("GenPVAssoc2RecoPV", "True signal vertex"),
    ("RecoAllAssoc2GenMatchedMerged", "Reconstructed vertices matched to a truth vertex, merged-vertex sample"),
    ("RecoAllAssoc2GenMatched", "Reconstructed vertices matched to a truth vertex"),
    ("RecoAllAssoc2GenMultiMatched", "Reconstructed vertices matched to more than one truth vertex (merged)"),
    ("RecoAllAssoc2MultiMatchedGen", "Reconstructed vertices whose truth vertex is also claimed by another reco vertex (split)"),
    ("RecoAllAssoc2GenSimForMerge", "Reconstructed vertices, merge-candidate truth pairing"),
    ("RecoAllAssoc2GenProperties", "Reconstructed vertex properties"),
    ("RecoAllAssoc2Gen", "All reconstructed vertices"),
    ("RecoAssoc2GenPVMatchedNotHighest", "Reconstructed vertices matched to the signal vertex but not ranked highest"),
    ("RecoAssoc2GenPVMatched", "Reconstructed vertices matched to the signal vertex"),
    ("RecoAssoc2GenPVNotMatched_GenPVTracksRemoved", "Reconstructed vertices not matched to the signal vertex, signal tracks removed"),
    ("RecoAssoc2GenPVNotMatched", "Reconstructed vertices not matched to the signal vertex"),
    ("RecoPVAssoc2GenPVMatched", "Reconstructed PV matched to the true signal vertex"),
    ("RecoPVAssoc2GenPVNotMatched", "Reconstructed PV not matched to the true signal vertex"),
]
PREFIX_LABELS.sort(key=lambda kv: -len(kv[0]))

EXPLICIT = {
    "TruePVLocationIndex": ("Signal-PV rank in the reconstructed vertex list", "Where the true signal vertex lands when reconstructed vertices are ranked by sum p_{T}^{2} (0 = correctly identified as the primary vertex)."),
    "TruePVLocationIndexCumulative": ("Signal-PV rank in the reconstructed vertex list (cumulative)", "Cumulative fraction of events where the true signal vertex ranks at or below a given position."),
    "TruePVLocationIndexSignalIsHighest": ("Signal-PV rank, events where signal has the highest p_{T}^{2}", "Where the true signal vertex lands in the ranked reco list, restricted to events where it genuinely has the highest sum p_{T}^{2}."),
    "TruePVLocationIndexSignalIsNotHighest": ("Signal-PV rank, events where signal does not have the highest p_{T}^{2}", "Where the true signal vertex lands in the ranked reco list, restricted to events where another vertex has higher sum p_{T}^{2}."),
    "KindOfSignalPV": ("Signal PV classification", "Categorizes how the true signal vertex was (mis)reconstructed: whether it matched any/the highest-ranked reco vertex."),
    "RecoAllAssoc2GenProperties": ("Reconstructed vertex classification", "Categorizes each reconstructed vertex by how it relates to truth vertices (matched, merged, split, fake, ...)."),
    "globalEfficiencies": ("Headline efficiency/fake-rate summary", "One bin per headline metric (overall efficiency, fake rate, etc.) integrated over the full sample."),
    "MisTagRate": ("Mis-tag rate", "Fraction of events where the true signal vertex was not ranked as the primary (first) reconstructed vertex."),
    "MisTagRateSignalIsHighest": ("Mis-tag rate, signal has highest p_{T}^{2}", "Mis-tag rate restricted to events where the signal vertex genuinely has the highest sum p_{T}^{2} of any truth vertex."),
    "MisTagRateSignalIsNotHighest": ("Mis-tag rate, signal is not highest p_{T}^{2}", "Mis-tag rate restricted to events where some other truth vertex has higher sum p_{T}^{2} than the signal."),
}

PULL_RESOL_RE = re.compile(r"^(.*_(Pull|Resol)[A-Za-z0-9]*)_vs_([A-Za-z0-9]+)_(Mean|Sigma)$")
RAW_PULL_RESOL_RE = re.compile(r"^(.*)_((?:Pull|Resol)[A-Za-z0-9]*)$")

STAT_DEFS = {
    "Pull": "(reconstructed minus true) divided by the estimated uncertainty; should look like a unit Gaussian if uncertainties are well calibrated.",
    "Resol": "Reconstructed minus true position; measures accuracy directly, in the same units as the quantity.",
}


def _prefix_context(name_or_quantity):
    for prefix, label in PREFIX_LABELS:
        if name_or_quantity.startswith(prefix + "_") or name_or_quantity == prefix:
            return label, name_or_quantity[len(prefix):].lstrip("_")
    return None, name_or_quantity


def _var_axis(var):
    return AXIS_LABELS.get(var, var)


def _var_word(var):
    return VAR_LABELS.get(var, var)


def derive(name):
    """Return (title, xlabel, ylabel, definition) for a DQM histogram name."""
    if name in EXPLICIT:
        title, definition = EXPLICIT[name]
        return title, "", "Events", definition

    for prefix, (ylab, definition) in RATE_DEFS.items():
        if name.startswith(prefix):
            var = name[len(prefix):]
            xlab = _var_axis(var)
            title = f"{ylab} vs. {_var_word(var)}"
            return title, xlab, ylab, definition

    if name.startswith("MisTagRate") and "_vs_" in name:
        base, var = name.split("_vs_", 1)
        base_title, base_def = EXPLICIT.get(base, ("Mis-tag rate", STAT_DEFS.get("", "")))
        xlab = _var_axis(var)
        title = f"{base_title} vs. {_var_word(var)}"
        return title, xlab, "Mis-tag rate", base_def

    if name.startswith("PV_purity_vs_"):
        var = name[len("PV_purity_vs_"):]
        xlab = _var_axis(var)
        title = f"PV finding efficiency vs. {_var_word(var)}"
        return title, xlab, "Efficiency", "Fraction of events where the true signal vertex has a matching reconstructed vertex."

    m = PULL_RESOL_RE.match(name)
    if m:
        quantity, kind, var, stat = m.group(1), m.group(2), m.group(3), m.group(4)
        statlabel = "Mean" if stat == "Mean" else "RMS"
        prefix_ctx, qsuffix = _prefix_context(quantity)
        qlabel = _var_word(qsuffix)
        xlab = _var_axis(var)
        ylab = f"{statlabel} of {qlabel}"
        ctx = f" ({prefix_ctx})" if prefix_ctx else ""
        title = f"{statlabel} of {qlabel} vs. {_var_word(var)}{ctx}"
        definition = STAT_DEFS.get(kind, "")
        return title, xlab, ylab, definition

    m = RAW_PULL_RESOL_RE.match(name)
    if m:
        base, qsuffix = m.group(1), m.group(2)
        prefix_ctx, _ = _prefix_context(base + "_")
        kind = "Pull" if qsuffix.startswith("Pull") else "Resol"
        qlabel = _var_word(qsuffix)
        title = f"{qlabel.capitalize()}" + (f" ({prefix_ctx})" if prefix_ctx else "")
        return title, AXIS_LABELS.get(qsuffix, qsuffix), "Vertices", STAT_DEFS.get(kind, "")

    for suffix in sorted(VAR_LABELS, key=len, reverse=True):
        if name.endswith("_" + suffix):
            base = name[: -(len(suffix) + 1)]
            prefix_ctx, _ = _prefix_context(base + "_")
            xlab = _var_axis(suffix)
            if prefix_ctx:
                title = f"{prefix_ctx}: {_var_word(suffix)}"
                definition = prefix_ctx[0].upper() + prefix_ctx[1:] + "."
            else:
                title = f"{_var_word(suffix).capitalize()} ({base})"
                definition = ""
            return title, xlab, "Vertices", definition

    return name.replace("_", " "), name, "Entries", ""
