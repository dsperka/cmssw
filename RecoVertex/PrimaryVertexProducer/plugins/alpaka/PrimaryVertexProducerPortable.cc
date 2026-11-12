#include "DataFormats/VertexSoA/interface/alpaka/VertexDeviceCollection.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexTrackDeviceCollection.h"
#include "DAInBlocksScratch.h"
#include "DataFormats/BeamSpot/interface/alpaka/BeamSpotDevice.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "TracksForDAInBlocksAlgo.h"
#include "DAInBlocksClusterizerAlgo.h"
#include "WeightedVertexFitterAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Reading nT() directly on host is undefined behavior on a real device (it's device
  // memory); only "works" by luck on SerialSync, where Device == DevHost.
  template <typename TColl>
  int32_t readNTOnHost(Queue& queue, TColl const& deviceTracks) {
    if constexpr (std::is_same_v<Device, alpaka_common::DevHost>) {
      return deviceTracks.view().nT();
    } else {
      auto hostTracks = cms::alpakatools::CopyToHost<TColl>::copyAsync(queue, deviceTracks);
      alpaka::wait(queue);
      return hostTracks.view().nT();
    }
  }

  /**
   * This class does vertexing by
   * - consuming set of Track
   * - clusterizing them into track clusters
   * - fitting cluster properties to vertex coordinates
   * - produces a device vertex product (Vertex)
   */
  class PrimaryVertexProducerPortable : public stream::EDProducer<> {
  public:
    PrimaryVertexProducerPortable(edm::ParameterSet const& config) : stream::EDProducer<>(config) {
      trackToken_ = consumes(config.getParameter<edm::InputTag>("TrackLabel"));
      beamSpotToken_ = consumes(config.getParameter<edm::InputTag>("BeamSpotLabel"));
      devicePutToken_ = produces();
      tracksInBlocksPutToken_ = produces();
      vertexTracksPutToken_ = produces();
      blockSize_ = config.getParameter<int32_t>("blockSize");
      blockOverlap_ = config.getParameter<double>("blockOverlap");
      // Fetch the (not free) ParameterSet once
      const auto tkClusParameters = config.getParameter<edm::ParameterSet>("TkClusParameters");
      clusterParams_ = {.Tmin = tkClusParameters.getParameter<double>("Tmin"),
                        .Tpurge = tkClusParameters.getParameter<double>("Tpurge"),
                        .Tstop = tkClusParameters.getParameter<double>("Tstop"),
                        .vertexSize = tkClusParameters.getParameter<double>("vertexSize"),
                        .coolingFactor = tkClusParameters.getParameter<double>("coolingFactor"),
                        .d0CutOff = tkClusParameters.getParameter<double>("d0CutOff"),
                        .dzCutOff = tkClusParameters.getParameter<double>("dzCutOff"),
                        .uniquetrkweight = tkClusParameters.getParameter<double>("uniquetrkweight"),
                        .uniquetrkminp = tkClusParameters.getParameter<double>("uniquetrkminp"),
                        .mintrkweight = tkClusParameters.getParameter<double>("mintrkweight"),
                        .zmerge = tkClusParameters.getParameter<double>("zmerge"),
                        .zrange = tkClusParameters.getParameter<double>("zrange"),
                        .convergence_mode = tkClusParameters.getParameter<int>("convergence_mode"),
                        .delta_lowT = tkClusParameters.getParameter<double>("delta_lowT"),
                        .delta_highT = tkClusParameters.getParameter<double>("delta_highT")};
      clusterParams_.uniquetrkminp = clusterParams_.uniquetrkminp * (1 - blockOverlap_);
      fitterParams_ = {
          .useBeamSpotConstraint =
              config.getParameter<edm::ParameterSet>("TkFitterParameters").getParameter<bool>("useBeamSpotConstraint"),
      };
    }

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
      const TrackForVertexDeviceCollection& inputtracks = iEvent.get(trackToken_);
      const BeamSpotDevice& beamSpot = iEvent.get(beamSpotToken_);
      int32_t nT = readNTOnHost(iEvent.queue(), inputtracks);
      int32_t nBlocks = nT > blockSize_ ? int32_t((nT - 1) / ((1. - blockOverlap_) * blockSize_))
                                        : 1;  // If all fit within a block, no need to split
      // Now the device collections we still need
      TrackForVertexDeviceCollection tracksInBlocks(iEvent.queue(), nBlocks * blockSize_);  // As high as needed
      TrackForVertexScratchCollection scratch(iEvent.queue(),
                                              nBlocks * blockSize_);  // device-only, not put in the event
      scratch.zeroInitialise(iEvent.queue());
      VertexScratchCollection vertexScratch(iEvent.queue(), 1024);  // device-only, same number of slots as deviceVertex
      vertexScratch.zeroInitialise(iEvent.queue());
      VertexDeviceCollection deviceVertex(
          iEvent.queue(), 1024);  // Hard capped to 1024, though we might want to restrict it for low PU cases

      VertexTrackDeviceCollection vertexTracks(iEvent.queue(),
                                               nBlocks * blockSize_);  // compact track lists of the vertices

      // run the algorithm
      //// First create the individual blocks
      TracksForDAInBlocksAlgo blockKernel_{};
      blockKernel_.createBlocks(iEvent.queue(), inputtracks, tracksInBlocks, blockSize_, blockOverlap_);
      // Need to have the blocks created before launching the next step
      //// Then run the clusterizer per blocks
      DAInBlocksClusterizerAlgo clusterizerKernel_{iEvent.queue(), blockSize_};
      clusterizerKernel_.clusterize(
          iEvent.queue(), tracksInBlocks, scratch, deviceVertex, vertexScratch, clusterParams_, nBlocks, blockSize_);
      clusterizerKernel_.resplit_tracks(
          iEvent.queue(), tracksInBlocks, scratch, deviceVertex, vertexScratch, clusterParams_, nBlocks, blockSize_);
      clusterizerKernel_.reject_outliers(
          iEvent.queue(), tracksInBlocks, scratch, deviceVertex, vertexScratch, clusterParams_, nBlocks, blockSize_);
      // Need to have all vertex before arbitrating and deciding what we keep
      clusterizerKernel_.arbitrate(
          iEvent.queue(), tracksInBlocks, scratch, deviceVertex, vertexScratch, clusterParams_, nBlocks, blockSize_);
      //// And then fit
      WeightedVertexFitterAlgo fitterKernel_{iEvent.queue(), fitterParams_};
      fitterKernel_.fit(iEvent.queue(), tracksInBlocks, deviceVertex, vertexScratch, vertexTracks, beamSpot);
      // Put the vertices in the event as a portable collection
      iEvent.emplace(devicePutToken_, std::move(deviceVertex));
      iEvent.emplace(tracksInBlocksPutToken_, std::move(tracksInBlocks));
      iEvent.emplace(vertexTracksPutToken_, std::move(vertexTracks));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("TrackLabel");
      desc.add<edm::InputTag>("BeamSpotLabel");
      desc.add<double>("blockOverlap");
      desc.add<int32_t>("blockSize");
      edm::ParameterSetDescription parf0;
      parf0.add<bool>("useBeamSpotConstraint", true);
      desc.add<edm::ParameterSetDescription>("TkFitterParameters", parf0);
      edm::ParameterSetDescription parc0;
      parc0.add<double>("d0CutOff", 3.0);
      parc0.add<double>("Tmin", 2.0);
      parc0.add<double>("delta_lowT", 0.001);
      parc0.add<double>("zmerge", 0.01);
      parc0.add<double>("dzCutOff", 3.0);
      parc0.add<double>("Tpurge", 2.0);
      parc0.add<int32_t>("convergence_mode", 0);
      parc0.add<double>("delta_highT", 0.01);
      parc0.add<double>("Tstop", 0.5);
      parc0.add<double>("coolingFactor", 0.6);
      parc0.add<double>("vertexSize", 0.006);
      parc0.add<double>("uniquetrkweight", 0.8);
      parc0.add<double>("uniquetrkminp", 0.0);
      parc0.add<double>("mintrkweight", 0.5);
      parc0.add<double>("zrange", 4.0);
      desc.add<edm::ParameterSetDescription>("TkClusParameters", parc0);
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    device::EDGetToken<TrackForVertexDeviceCollection> trackToken_;
    device::EDGetToken<BeamSpotDevice> beamSpotToken_;
    device::EDPutToken<VertexDeviceCollection> devicePutToken_;
    device::EDPutToken<TrackForVertexDeviceCollection> tracksInBlocksPutToken_;
    device::EDPutToken<VertexTrackDeviceCollection> vertexTracksPutToken_;
    int32_t blockSize_;
    double blockOverlap_;
    FitterParameters fitterParams_;
    DAInBlocksClusterParameters clusterParams_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(PrimaryVertexProducerPortable);
