"""Human-readable titles, axis labels and descriptions for the histograms of the MTD vertex validation
(Validation/MtdValidation, Primary4DVertexValidation: DQM folder MTD/Vertices/<collection>/).

derive(name, own) -> (title, xlabel, ylabel, definition), or None when the name is not recognised.
  title       plain text for the gallery card
  xlabel      ROOT (TLatex) axis title, kept short enough for the plot
  ylabel      ROOT axis title
  definition  one or two sentences for the gallery card
`own` is an optional dict with the histogram's own title / axis titles from the DQM file.

What the validation does (read from Primary4DVertexValidation.cc), as used in the descriptions below:
 * Every reconstructed vertex is matched to a simulated vertex through the tracks it contains. Only tracks with a vertex weight of at
   least `trackweightTh` (0.5) count; tracks need |eta| <= 3 and pT > 0.7 GeV. A track's weight "wos" is its vertex weight divided by its z
   variance (dz^2 plus the beam-spot contribution plus a 20 um floor), and further divided by erf(sigma_t / sigma_T) for tracks with MTD time
   information; "wnt" is the vertex weight times min(pT, 1 GeV).
 * A reconstructed vertex is "dominated" by the simulated vertex with the largest summed wos (at least two tracks). The one-to-one matching
   runs in ranks 1..8: rank 1 for simulated vertices that dominate a single reconstructed vertex, higher ranks for ambiguous cases.
   A reconstructed vertex with rank 1..7 is "real"; one without a match (or only the weakest rank) is "fake". A fake vertex that still
   takes more than 30% of its wos from one simulated vertex is a "split" (that vertex was reconstructed more than once); the others are
   "other fake".
 * Only vertices with ndof > 4 are classified and plotted.
 * "LV" / "RecLV" is the leading reconstructed vertex (the one matched to the simulated signal vertex); histograms without it use all
   matched reconstructed vertices.
 * "PU tracks" come from another simulated vertex than the matched one, "fake tracks" have no associated simulated particle,
   "secondary tracks" come from secondary vertices.
"""
import re

TRK_CLASS = {"PUTrack": "pileup", "FakeTrack": "fake", "SecTrack": "secondary-vertex"}
QTY = {
    "Mult": ("number of tracks", "number of tracks"),
    "SumPt": ("sum of track pT", "sum of track p_{T}"),
    "SumPt2": ("sum of track pT squared", "sum of track p_{T}^{2}"),
    "SumWnt": ("weighted number of tracks (weight x min(pT, 1 GeV))", "weighted number of tracks"),
    "SumWos": ("sum of track weights over z variance (wos)", "sum of track weights over variance"),
    "Wnt": ("weighted track count (weight x min(pT, 1 GeV))", "weighted track count"),
}
JET_QTY = {
    "Mult": ("jet multiplicity", "number of jets"),
    "Ht": ("scalar sum of jet ET (HT)", "HT"),
    "MetPt": ("missing transverse energy (MET)", "MET"),
    "SumPt2": ("sum of jet pT squared", "sum of jet p_{T}^{2}"),
    "SumPz": ("sum of jet pz", "sum of jet p_{z}"),
}
SPECIES = {"K": "kaon", "P": "proton", "Pi": "pion"}
REGION = {"Barrel": ("barrel timing layer", "|#eta| < 1.5"), "Endcap": ("endcap timing layer", "1.6 < |#eta| < 3")}
ETA_PART = {"_lowEta": "1.6 < |#eta| < 2.1", "_highEta": "2.1 < |#eta| < 3"}
MVA_CAT = {"LowMVA": "time-quality MVA < 0.5", "MediumMVA": "time-quality MVA between 0.5 and 0.8", "HighMVA": "time-quality MVA > 0.8"}
MOM_CAT = {"HighP": "momentum above 2 GeV", "LowP": "momentum below 2 GeV"}


def _plain(s):
    return s.replace("|#eta|", "|eta|").replace("#eta", "eta")


def _pid(name):
    m = re.match(r"^(Barrel|Endcap)(.+?)(_lowEta|_highEta)?$", name)
    if not m:
        return None
    reg, body, part = m.group(1), m.group(2), m.group(3)
    rname, reta = REGION[reg]
    eta = ETA_PART.get(part, reta)
    where = f"{rname}, {eta}"
    xp = "Track momentum p [GeV]"
    for pat, builder in (
        (r"True(K|P|Pi)As(K|P|Pi)", lambda a: (f"True {SPECIES[a[0]]}s identified as {SPECIES[a[1]]}s: momentum spectrum",
                                              xp, "Tracks", f"Momentum of tracks of the leading vertex that are truly {SPECIES[a[0]]}s and were assigned the {SPECIES[a[1]]} mass hypothesis.")),
        (r"True(K|P|Pi)NoPID", lambda a: (f"True {SPECIES[a[0]]}s with no PID decision: momentum spectrum", xp, "Tracks",
                                          f"Momentum of tracks that are truly {SPECIES[a[0]]}s but received no particle identification (see the PID failure categories).")),
        (r"True(K|P|Pi)", lambda a: (f"True {SPECIES[a[0]]}s: momentum spectrum", xp, "Tracks",
                                     f"Momentum of tracks of the leading vertex that are truly {SPECIES[a[0]]}s (denominator of the identification fractions).")),
        (r"As(K|P|Pi)", lambda a: (f"Tracks identified as {SPECIES[a[0]]}s: momentum spectrum", xp, "Tracks",
                                   f"Momentum of tracks assigned the {SPECIES[a[0]]} mass hypothesis, whatever their true species.")),
        (r"NoPIDtype", lambda a: ("Why no PID decision was made", "PID failure category", "Tracks",
                                  "1: no PID probabilities available; 2: probability not finite; 3: fall-back to the pion hypothesis because the PID-based time uncertainty is smaller than the safe one.")),
        (r"PIDp", lambda a: ("Tracks with a PID decision: momentum spectrum", xp, "Tracks", "Momentum of leading-vertex tracks that have MTD time and received a particle identification.")),
        (r"NoPID", lambda a: ("Tracks with no PID decision: momentum spectrum", xp, "Tracks", "Momentum of tracks that received no particle identification.")),
        (r"PID(K|P|Pi)NoPIDEff", lambda a: (f"Fraction of true {SPECIES[a[0]]}s with no PID decision vs momentum", xp, "Fraction of true tracks",
                                            f"Fraction of true {SPECIES[a[0]]}s that received no particle identification, as a function of momentum.")),
        (r"PID(K|P|Pi)As(K|P|Pi)Eff", lambda a: (f"Fraction of true {SPECIES[a[0]]}s identified as {SPECIES[a[1]]}s vs momentum", xp,
                                                 "Fraction of true tracks" if a[0] != a[1] else "Identification efficiency",
                                                 f"Fraction of true {SPECIES[a[0]]}s assigned the {SPECIES[a[1]]} hypothesis, as a function of momentum." + (" This is the identification efficiency." if a[0] == a[1] else " This is a misidentification rate."))),
        (r"(K|P|Pi)Purity", lambda a: (f"Purity of {SPECIES[a[0]]} identification vs momentum", xp, f"Fraction truly {SPECIES[a[0]]}s",
                                       f"Among the tracks identified as {SPECIES[a[0]]}s, the fraction that truly are {SPECIES[a[0]]}s, as a function of momentum.")),
    ):
        mm = re.fullmatch(pat, body)
        if mm:
            title, xl, yl, d = builder(mm.groups())
            return f"{title} ({_plain(where)})", xl, yl, d + f" Tracks in the {_plain(where)}."
    return None


def _track_composition(name):
    m = re.match(r"^(PUTrack|FakeTrack|SecTrack)(RecLV)?(Rel)?(Mult|SumPt2|SumPt|SumWnt|SumWos|Wnt)(vs(Mult|SumPt2|SumPt|SumWnt|SumWos))?$", name)
    if not m:
        return None
    cls, reclv, rel, qty, vs = TRK_CLASS[m.group(1)], m.group(2), m.group(3), m.group(4), m.group(5)
    qplain, qroot = QTY[qty]
    scope = "the leading reconstructed vertex" if reclv else "all matched reconstructed vertices"
    scope_short = "leading vertex" if reclv else "matched vertex"
    ev = "Events" if reclv else "Matched vertices"
    if vs:
        vq = m.group(6)
        xq = QTY[vq][1]
        xl = "Number of tracks in the vertex" if vq == "Mult" else f"log10({xq} of the vertex tracks)"
        return (f"Fraction of {cls} tracks vs size of the {scope_short} ({qplain})", xl,
                f"Mean fraction of {cls} tracks" if len(cls) < 10 else "Mean fraction of tracks",
                f"Average fraction of {cls} tracks (by {qplain}) in {scope}, as a function of the vertex size.")
    if rel:
        xl = f"Fraction of {cls} tracks in the {scope_short}"
        if qty != "Mult":
            xl += f" (by {qroot})"
        return (f"Fraction of {cls} tracks in {scope}, by {qplain}", xl, ev,
                f"Share of the {qplain} of the vertex that comes from {cls} tracks (0 = none, 1 = all). One entry per {'event' if reclv else 'matched vertex'}.")
    if qty == "Mult":
        xl = f"Number of {cls} tracks in the {scope_short}"
    elif qty == "Wnt" or True:
        xl = f"log10({cls} tracks: {qroot})" if qty != "Wnt" else f"{cls.capitalize()} track weight x min(pT, 1 GeV)"
    return (f"{qplain[0].upper() + qplain[1:]} of {cls} tracks in {scope}", xl, ev,
            f"{qplain[0].upper() + qplain[1:]} of the {cls} tracks in {scope}.")


def _jets(name):
    m = re.match(r"^Jets(RecLV)?(PU|Fake)(Rel)?(Mult|Ht|MetPt|SumPt2|SumPz)(vs(Mult|Ht|MetPt|SumPt2|SumPz))?$", name)
    if not m:
        return None
    reclv, kind, rel, qty, vs = m.group(1), m.group(2), m.group(3), m.group(4), m.group(5)
    qplain, qroot = JET_QTY[qty]
    scope = "leading vertex" if reclv else "matched vertices"
    what = "pileup tracks" if kind == "PU" else "fake tracks"
    ev = "Events" if reclv else "Matched vertices"
    if vs:
        return (f"Relative {what} contribution to {qplain} vs its size ({scope})", f"log10({JET_QTY[m.group(6)][1]})" if m.group(6) != "Mult" else "Number of jets",
                "Mean relative contribution", f"Average relative contribution of {what} to the {qplain} of jets built from the tracks of the {scope}.")
    if rel:
        return (f"Relative contribution of {what} to {qplain} ({scope})", f"Relative contribution of {what} to {qroot}", ev,
                f"(value with all tracks - value without {what}) / value with all tracks, for jets built from the tracks of the {scope}.")
    return (f"Contribution of {what} to {qplain} ({scope})", f"log10(contribution of {what} to {qroot})" if qty != "Mult" else f"Jets lost without {what}", ev,
            f"Difference of the {qplain} of jets built from the tracks of the {scope} with and without the {what}.")


def _simple(name):
    t = {
        "RecVerNumber": ("Reconstructed vertices per event", "Reconstructed vertices in the event", "Events", "All reconstructed vertices of the collection in an event, before any quality selection."),
        "RecSelVerNumber": ("Reconstructed vertices per event passing ndof > 4", "Reconstructed vertices with ndof > 4 (real + fake)", "Events",
                            "Vertices of the collection with more than 4 degrees of freedom, the ones that are classified as real or fake."),
        "SimVerNumber": ("Simulated vertices per event", "Simulated vertices in the event (pileup + signal)", "Events", "Number of simulated primary vertices in the event."),
        "MatchQuality": ("Quality of the reconstructed-to-simulated vertex match", "Match rank (1 = unambiguous, up to 8 = weakest)", "Reconstructed vertices",
                         "Rank at which a reconstructed vertex was matched one-to-one to a simulated vertex. Lower ranks are cleaner matches."),
        "PUvsReal": ("Real vertices per event vs pileup", "Simulated vertices in the event (pileup)", "Mean real reconstructed vertices", "Average number of reconstructed vertices classified as real (matched to a simulated vertex) per event."),
        "PUvsFake": ("Fake vertices per event vs pileup", "Simulated vertices in the event (pileup)", "Mean fake reconstructed vertices", "Average number of reconstructed vertices without a simulated match per event (includes splits)."),
        "PUvsSplit": ("Split vertices per event vs pileup", "Simulated vertices in the event (pileup)", "Mean split vertices", "Average number of fake vertices that are copies of a real one (they take over 30% of their weight from a single simulated vertex)."),
        "PUvsOtherFake": ("Other fake vertices per event vs pileup", "Simulated vertices in the event (pileup)", "Mean other fake vertices", "Average number of fake vertices that are not splits of a simulated vertex."),
        "RecoPVPosSignal": ("Rank of the signal vertex in the reconstructed collection", "Position of the signal-matched vertex in the collection", "Events", "Where, in the collection sorted by sum pT squared, the vertex matched to the simulated signal sits (0 = leading)."),
        "RecoPVPosSignalNotHighestPt": ("Rank of the signal vertex when it is not the leading one", "Position of the signal-matched vertex in the collection", "Events", "Same as above, only for events where the signal vertex is not the first."),
        "RecoPosInRecoOrigCollection": ("Index of the signal vertex in the original reconstructed collection", "Reconstructed vertex index", "Events", "Index of the reconstructed signal vertex in the collection before sorting."),
        "RecoPosInSimCollection": ("Index of the simulated vertex matched to the reconstructed signal", "Simulated vertex index", "Events", "Index in the simulated collection of the vertex matched to the reconstructed signal vertex."),
        "SimPosInSimOrigCollection": ("Index of the simulated signal vertex", "Simulated vertex index", "Events", "Index of the simulated signal vertex in the original simulated collection."),
        "recPVZ": ("z of reconstructed vertices", "Reconstructed vertex z [cm]", "Vertices per 10 mm", "z position of the reconstructed vertices with ndof > 4."),
        "recPVT": ("Time of reconstructed vertices", "Reconstructed vertex time [ns]", "Vertices per 50 ps", "Time of the reconstructed vertices that have a time measurement."),
        "simPVZ": ("z of simulated vertices", "Simulated vertex z [cm]", "Vertices per 10 mm", "z position of the simulated vertices."),
        "simPVT": ("Time of simulated vertices", "Simulated vertex time [ns]", "Vertices per 50 ps", "Time of the simulated vertices."),
        "simPVTvsZ": ("Time of simulated vertices vs z", "Simulated vertex z [cm]", "Mean simulated vertex time [ns]", "Average simulated vertex time as a function of z (collision time spread along the beam line)."),
        "TimeRes": ("Vertex time residual", "Reconstructed minus simulated vertex time [ns]", "Matched vertices", "Reconstructed vertex time minus the time of its matched simulated vertex, for matched vertices."),
        "TimeSignalRes": ("Signal vertex time residual", "Reconstructed minus simulated signal vertex time [ns]", "Events", "Same as the vertex time residual, for the signal vertex only."),
        "TimePull": ("Vertex time pull", "Time residual / reconstructed time uncertainty", "Matched vertices", "Time residual divided by the uncertainty of the reconstructed vertex time. A unit-width Gaussian means correct uncertainties."),
        "TimeSignalPull": ("Signal vertex time pull", "Signal time residual / reconstructed time uncertainty", "Events", "Time pull of the signal vertex only."),
        "SelVtxTrk": ("Tracks selected for vertexing", "log10(tracks selected for the primary vertex)", "Events", "Number of tracks passing the vertexing track selection, per event."),
        "SelVtxTrkwTime": ("Selected tracks that have MTD time", "log10(selected tracks with MTD time)", "Events", "Number of vertexing tracks with a time measurement from the MTD, per event."),
        "SelVtxTrkvsEta": ("Pseudorapidity of the selected tracks", "Track pseudorapidity #eta", "Tracks", "Distribution in eta of the tracks selected for vertexing."),
        "SelVtxTrkwTimevsEta": ("Pseudorapidity of the selected tracks with MTD time", "Track pseudorapidity #eta", "Tracks", "Distribution in eta of the selected tracks that have MTD time."),
        "UnAssocTracks": ("Tracks not associated to a simulated particle", "log10(unassociated tracks in the event)", "Events", "Number of reconstructed tracks without an associated simulated particle."),
        "UnAssocTracksFake": ("Unassociated tracks that are fake", "log10(unassociated fake tracks in the event)", "Events", "Number of tracks without an associated simulated particle that are classified as fake."),
        "FractionUnAssocTracks": ("Fraction of unassociated tracks", "Fraction of tracks without a simulated match", "Events", "Fraction of the reconstructed tracks of the event that have no associated simulated particle."),
        "FractionUnAssocTracksFake": ("Fraction of unassociated fake tracks", "Fraction of unassociated tracks that are fake", "Events", "Fraction of the tracks without a simulated match that are classified as fake."),
        "EffPtTot": ("pT of the leading-vertex tracks", "Track p_{T} [GeV]", "Tracks", "Transverse momentum of the tracks associated to the leading vertex (denominator of the time-measurement efficiency)."),
        "EffEtaTot": ("Pseudorapidity of the leading-vertex tracks", "Track |#eta|", "Tracks", "Absolute pseudorapidity of the tracks associated to the leading vertex."),
        "MatchedTPEffPtTot": ("pT of leading-vertex tracks matched to a simulated particle", "Track p_{T} [GeV]", "Tracks", "pT of the leading-vertex tracks that are matched to a simulated particle."),
        "MatchedTPEffPtMtd": ("pT of matched leading-vertex tracks with MTD time", "Track p_{T} [GeV]", "Tracks", "Same, only tracks that also have a time measurement."),
        "MatchedTPEffEtaTot": ("Pseudorapidity of leading-vertex tracks matched to a simulated particle", "Track |#eta|", "Tracks", "|eta| of the leading-vertex tracks that are matched to a simulated particle."),
        "MatchedTPEffEtaMtd": ("Pseudorapidity of matched leading-vertex tracks with MTD time", "Track |#eta|", "Tracks", "Same, only tracks that also have a time measurement."),
        "TPPtMatchEff": ("Efficiency to match tracks to simulated particles vs pT", "Particle p_{T} [GeV]", "Efficiency", "Fraction of the selected simulated particles of the leading vertex that are matched to a reconstructed track, vs pT."),
        "TPEtaMatchEff": ("Efficiency to match tracks to simulated particles vs eta", "Particle |#eta|", "Efficiency", "Same, vs pseudorapidity."),
        "TPPtSelEff": ("Efficiency of the track selection vs pT", "Particle p_{T} [GeV]", "Efficiency", "Fraction of the leading-vertex simulated particles whose track passes the selection, vs pT."),
        "TPEtaSelEff": ("Efficiency of the track selection vs eta", "Particle |#eta|", "Efficiency", "Same, vs pseudorapidity."),
        "MatchedTPTrackRes": ("Track time residual (leading vertex, matched)", "Track time minus simulated time [ns]", "Tracks", "Reconstructed track time minus the time of the matched simulated particle, for matched tracks of the leading vertex."),
        "MatchedTPTrackPull": ("Track time pull (leading vertex, matched)", "Track time residual / its uncertainty", "Tracks", "Time residual divided by the track time uncertainty."),
        "MatchedTPTrackZposResTot": ("Track z residual (leading vertex, matched)", "Track z at closest approach minus simulated z [cm]", "Tracks", "Difference in z between the track at its point of closest approach and the simulated production point."),
        "TrackZposResTot": ("Track z residual", "Track z at closest approach minus simulated z [cm]", "Tracks", "Difference in z between the track at its point of closest approach and the simulated production point."),
        "VtxTrackMult": ("Tracks per vertex", "log10(tracks in the vertex)", "Vertices", "Number of tracks in each reconstructed vertex."),
        "VtxTrackMultLin": ("Tracks per vertex", "Tracks in the vertex", "Vertices", "Number of tracks in each reconstructed vertex."),
        "VtxTrackMultPassNdof": ("Tracks per vertex, ndof > 4", "log10(tracks in the vertex)", "Vertices", "Track count of vertices with more than 4 degrees of freedom."),
        "VtxTrackMultPassNdofLin": ("Tracks per vertex, ndof > 4", "Tracks in the vertex", "Vertices", "Track count of vertices with more than 4 degrees of freedom."),
        "VtxTrackMultFailNdof": ("Tracks per vertex, ndof <= 4", "Tracks in the vertex", "Vertices", "Track count of vertices that fail the ndof > 4 requirement."),
        "VtxTrackRecLVMult": ("Tracks in the leading vertex", "log10(tracks in the leading vertex)", "Events", "Number of tracks of the leading reconstructed vertex (those with vertex weight above 0.5)."),
        "VtxTrackRecLVMultLin": ("Tracks in the leading vertex", "Tracks in the leading vertex", "Events", "Number of tracks of the leading reconstructed vertex (those with vertex weight above 0.5)."),
        "VtxTrackW": ("Track weights in vertices", "Track weight in the vertex", "Tracks", "Weight each track has in the vertex it belongs to (1 = fully assigned)."),
        "VtxTrackRecLVW": ("Track weights in the leading vertex", "Track weight in the leading vertex", "Tracks", "Weight each track has in the leading vertex."),
        "VtxTrackWnt": ("Weighted track count in vertices", "Track weight x min(pT, 1 GeV)", "Tracks", "Track weight times min(pT, 1 GeV), per track, for all vertices."),
        "VtxTrackRecLVWnt": ("Weighted track count in the leading vertex", "Track weight x min(pT, 1 GeV)", "Tracks", "Track weight times min(pT, 1 GeV), per track, for the leading vertex."),
    }
    return t.get(name)


def _delta(name):
    m = re.match(r"^Delta(Z|T)(realreal|fakefake|fakereal)$", name)
    if not m:
        return None
    var, pair = m.group(1), m.group(2)
    kinds = {"realreal": "two real vertices", "fakefake": "two fake vertices", "fakereal": "a fake and a real vertex"}[pair]
    if var == "Z":
        return (f"Distance in z between {kinds}", "|z difference| between the two vertices [cm]", "Vertex pairs", f"|#Delta z| for pairs made of {kinds} in the same event. Close pairs of fake and real vertices indicate splitting.")
    return (f"Distance in time between {kinds}", "|time difference| / its uncertainty", "Vertex pairs", f"|#Delta t| in units of its uncertainty for pairs made of {kinds} in the same event.")


def _tracks_time(name):
    m = re.match(r"^Track(Res|Pull)(Mass(?:True)?)?(HighP|LowP)?(?:-(Pions|Protons))?(?:-(Low|Medium|High)MVA)?$", name)
    if m and not (m.group(2) and m.group(3)):
        kind, mass, mom, spec, mva = m.group(1), m.group(2), m.group(3), m.group(4), m.group(5)
        subset = []
        if spec:
            subset.append(spec.lower()[:-1] + " tracks" if spec == "Pions" else "proton tracks")
        if mom:
            subset.append(MOM_CAT[mom])
        if mva:
            subset.append(MVA_CAT[mva + "MVA"])
        what = ("tracks" if not spec else f"{spec.lower()[:-1]} tracks") + (f" with {', '.join(s for s in subset if 'tracks' not in s)}" if any('tracks' not in s for s in subset) else "")
        if kind == "Res" and not mass:
            return (f"Track time residual: {what}", "Track time minus simulated time [ns]", "Tracks", f"Reconstructed MTD track time minus the simulated time, for {what} of the leading vertex.")
        if kind == "Res" and mass == "Mass":
            return (f"Track time vs mass-hypothesis estimate: {what}", "Track time minus estimate from the true mass [ns]", "Tracks",
                    f"Reconstructed track time minus the time estimated from the MTD hit time and the time of flight under the true particle mass, for {what}.")
        if kind == "Res" and mass == "MassTrue":
            return (f"Mass-hypothesis time estimate vs simulation: {what}", "Estimate from the true mass minus simulated time [ns]", "Tracks",
                    f"Time estimated from the MTD hit time and the time of flight under the true mass, minus the simulated time, for {what}. Isolates the effect of the time-of-flight correction.")
        return (f"Track time pull: {what}", "Track time residual / its uncertainty", "Tracks", f"Track time residual divided by its uncertainty for {what}. A unit-width Gaussian means correct uncertainties.")
    m = re.match(r"^Track(Z|3D)posRes(Tot)?(?:-(Low|Medium|High)MVA)?$", name)
    if m:
        dim, tot, mva = m.group(1), m.group(2), m.group(3)
        what = "tracks" + (f" with {MVA_CAT[mva + 'MVA']}" if mva else "")
        if dim == "Z":
            return (f"Track z residual: {what}", "Track z at closest approach minus simulated z [cm]", "Tracks", f"z difference between the track and the simulated particle for {what}.")
        return (f"Track 3D position residual: {what}", "3D distance of track to simulated point [cm]", "Tracks", f"3D distance between the track's closest approach and the simulated particle origin, for {what}.")
    m = re.match(r"^(NoPIDTrack)?Track(MVA|TimeRes|TimePull|TimeSigma)(NoPID|CorrectPID|WrongPID)(?:true(K|P|Pi)|Type([123]))?$", name) or re.match(r"^(NoPIDTrack)(MVA|Sigma|TimePull|TimeRes)NoPIDType([123])$", name)
    if m:
        g = m.groups()
        q = [x for x in g if x in ("MVA", "TimeRes", "TimePull", "TimeSigma", "Sigma")][0]
        if g[0] == "NoPIDTrack" and len(g) == 3:
            pidtxt, sub = "no PID", f"PID failure category {g[2]}"
        else:
            pidtxt = {"NoPID": "no PID", "CorrectPID": "correctly identified", "WrongPID": "wrongly identified"}[g[2]]
            sub = (f"true {SPECIES[g[3]]}s" if g[3] else (f"PID failure category {g[4]}" if g[4] else ""))
        what = f"tracks with {pidtxt}" + (f" ({sub})" if sub else "")
        if q == "MVA":
            return (f"Time-quality MVA of {what}", "MTD track time-quality MVA score", "Tracks", f"Score of the MVA that rates the quality of the MTD time of {what}.")
        if q in ("TimeSigma", "Sigma"):
            return (f"Track time uncertainty of {what}", "Track time uncertainty [ns]", "Tracks", f"Uncertainty of the track time of {what}.")
        if q == "TimeRes":
            return (f"Track time residual of {what}", "Track time minus simulated time [ns]", "Tracks", f"Reconstructed track time minus simulated time for {what}.")
        return (f"Track time pull of {what}", "Track time residual / its uncertainty", "Tracks", f"Time residual divided by its uncertainty for {what}.")
    return None


def derive(name, own=None):
    for fn in (_pid, _simple, _track_composition, _jets, _delta, _tracks_time):
        r = fn(name)
        if r:
            return r
    return None
