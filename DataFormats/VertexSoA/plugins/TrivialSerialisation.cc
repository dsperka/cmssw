#include "DataFormats/VertexSoA/interface/ZVertexHost.h"
#include "DataFormats/VertexSoA/interface/VertexHostCollection.h"
#include "DataFormats/VertexSoA/interface/TrackForVertexHostCollection.h"
#include "DataFormats/VertexSoA/interface/VertexTrackHostCollection.h"
#include "HeterogeneousCore/TrivialSerialisation/interface/SerialiserFactory.h"

DEFINE_TRIVIAL_SERIALISER_PLUGIN(reco::ZVertexHost);
DEFINE_TRIVIAL_SERIALISER_PLUGIN(::reco::VertexHostCollection);
DEFINE_TRIVIAL_SERIALISER_PLUGIN(::reco::TrackForVertexHostCollection);
DEFINE_TRIVIAL_SERIALISER_PLUGIN(::reco::VertexTrackHostCollection);
