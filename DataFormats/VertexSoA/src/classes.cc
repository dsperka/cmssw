#include "DataFormats/Portable/interface/PortableHostCollectionReadRules.h"
#include "DataFormats/VertexSoA/interface/ZVertexHost.h"
#include "DataFormats/VertexSoA/interface/VertexGNNHostCollection.h"
#include "DataFormats/VertexSoA/interface/VertexHostCollection.h"
#include "DataFormats/VertexSoA/interface/TrackForVertexHostCollection.h"
#include "DataFormats/VertexSoA/interface/VertexTrackHostCollection.h"

SET_PORTABLEHOSTCOLLECTION_READ_RULES(reco::ZVertexHost);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(::reco::VertexHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(::reco::TrackForVertexHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(vertexgnn::TrackFeaturesHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(vertexgnn::GNNOutputHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(vertexgnn::GNNSlotHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(::reco::VertexTrackHostCollection);
