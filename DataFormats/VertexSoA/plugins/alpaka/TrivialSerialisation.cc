#include "DataFormats/VertexSoA/interface/VertexGNNHostCollection.h"
#include "DataFormats/VertexSoA/interface/ZVertexHost.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexGNNDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/TrackForVertexHostCollection.h"
#include "DataFormats/VertexSoA/interface/VertexHostCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/TrackForVertexDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/VertexTrackHostCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexTrackDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/ZVertexSoACollection.h"
#include "HeterogeneousCore/TrivialSerialisation/interface/alpaka/SerialiserFactoryDevice.h"

DEFINE_TRIVIAL_SERIALISER_PORTABLE_PLUGIN(reco::ZVertexHost, reco::ZVertexSoACollection);
DEFINE_TRIVIAL_SERIALISER_PORTABLE_PLUGIN(vertexgnn::TrackFeaturesHostCollection,
                                          vertexgnn::TrackFeaturesDeviceCollection);
DEFINE_TRIVIAL_SERIALISER_PORTABLE_PLUGIN(vertexgnn::GNNOutputHostCollection, vertexgnn::GNNOutputDeviceCollection);
DEFINE_TRIVIAL_SERIALISER_PORTABLE_PLUGIN(vertexgnn::GNNSlotHostCollection, vertexgnn::GNNSlotDeviceCollection);
DEFINE_TRIVIAL_SERIALISER_PORTABLE_PLUGIN(::reco::VertexHostCollection, VertexDeviceCollection);
DEFINE_TRIVIAL_SERIALISER_PORTABLE_PLUGIN(::reco::TrackForVertexHostCollection, TrackForVertexDeviceCollection);
DEFINE_TRIVIAL_SERIALISER_PORTABLE_PLUGIN(::reco::VertexTrackHostCollection, VertexTrackDeviceCollection);
