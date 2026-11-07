"""Experiment: merge two adjacent vertices only if closer than min(mergeSignificance*combinedSigma, 2*(2*vertexSize)).
Old rule: max(2*vertexSize, 3*combinedSigma) -- the 0.012 cm floor over-merges at high track multiplicity."""
import pathlib
R = pathlib.Path("/home/claude/scratch/ws2/CMSSW_16_1_0_pre4/src/RecoVertex/PrimaryVertexProducer")
p = R/"src/DAClusterizerInZ_vect.cc"; s = p.read_text()
a = """    double mergeThreshold = 2 * vertexSize_;
    if (k + 1 < nv) {
      double combinedSigma = std::sqrt(vtx_sigma_z[k] * vtx_sigma_z[k] + vtx_sigma_z[k + 1] * vtx_sigma_z[k + 1]);
      mergeThreshold = std::max(mergeThreshold, mergeSignificance * combinedSigma);
    }"""
b = """    double mergeThreshold = 2 * vertexSize_;
    if (k + 1 < nv) {
      double combinedSigma = std::sqrt(vtx_sigma_z[k] * vtx_sigma_z[k] + vtx_sigma_z[k + 1] * vtx_sigma_z[k + 1]);
      // z resolution improves with the number of tracks: use the statistical significance, but never merge
      // vertices that are further apart than 2 * (2 * vertexSize)
      mergeThreshold = std::min(mergeSignificance * combinedSigma, 2. * mergeThreshold);
    }"""
assert s.count(a) == 1; p.write_text(s.replace(a, b))
q = R/"plugins/alpaka/DAInBlocksClusterizerArbitrator.dev.cc"; t = q.read_text()
c = "double mergeThreshold = std::max(2 * cParams.vertexSize, 3.0 * combinedSigma);"
d = "double mergeThreshold = std::min(3.0 * combinedSigma, 4 * cParams.vertexSize);"
assert t.count(c) == 1; q.write_text(t.replace(c, d)); print("merge patch applied")
