// -*- C++ -*-
//
// Package:    HLTrigger/JetMET
// Class:      HLTScoutingPFProducer
//
/**\class HLTScoutingPFProducer HLTScoutingPFProducer.cc HLTrigger/JetMET/plugins/HLTScoutingPFProducer.cc

Description: Producer for ScoutingPFJets from reco::PFJet objects, ScoutingVertexs from reco::Vertexs and ScoutingParticles from reco::PFCandidates

*/
//
// Original Author:  Dustin James Anderson
//         Created:  Fri, 12 Jun 2015 15:49:20 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticleV2.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Math/interface/libminifloat.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/BTauReco/interface/DeepBoostedJetTagInfo.h"

using namespace btagbtvdeep;

class HLTScoutingPFProducerV2 : public edm::global::EDProducer<> {
public:
  explicit HLTScoutingPFProducerV2(const edm::ParameterSet &);
  ~HLTScoutingPFProducerV2() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::StreamID sid, edm::Event &iEvent, edm::EventSetup const &setup) const final;

  const edm::EDGetTokenT<reco::PFJetCollection> pfJetCollection_;
  const edm::EDGetTokenT<reco::JetTagCollection> pfJetTagCollection_;
  const edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidateCollection_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
  const edm::EDGetTokenT<reco::PFMETCollection> metCollection_;
  const edm::EDGetTokenT<double> rho_;

  const double pfJetPtCut;
  const double pfJetEtaCut;
  const double pfCandidatePtCut;
  const double pfCandidateEtaCut;
  const int mantissaPrecision;

  const bool doJetTags;
  const bool doCandidates;
  const bool doMet;

  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> trackBuilderToken;
};

//
// constructors and destructor
//
HLTScoutingPFProducerV2::HLTScoutingPFProducerV2(const edm::ParameterSet &iConfig)
    : pfJetCollection_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("pfJetCollection"))),
      pfJetTagCollection_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("pfJetTagCollection"))),
      pfCandidateCollection_(
          consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidateCollection"))),
      vertexCollection_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
      metCollection_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metCollection"))),
      rho_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
      pfJetPtCut(iConfig.getParameter<double>("pfJetPtCut")),
      pfJetEtaCut(iConfig.getParameter<double>("pfJetEtaCut")),
      pfCandidatePtCut(iConfig.getParameter<double>("pfCandidatePtCut")),
      pfCandidateEtaCut(iConfig.getParameter<double>("pfCandidateEtaCut")),
      mantissaPrecision(iConfig.getParameter<int>("mantissaPrecision")),
      doJetTags(iConfig.getParameter<bool>("doJetTags")),
      doCandidates(iConfig.getParameter<bool>("doCandidates")),
      doMet(iConfig.getParameter<bool>("doMet")),
      trackBuilderToken(
        esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))) {
  //register products
  produces<Run3ScoutingPFJetCollection>();
  produces<Run3ScoutingParticleV2Collection>();
  produces<double>("rho");
  produces<double>("pfMetPt");
  produces<double>("pfMetPhi");
}

HLTScoutingPFProducerV2::~HLTScoutingPFProducerV2() = default;

// ------------ method called to produce the data  ------------
void HLTScoutingPFProducerV2::produce(edm::StreamID sid, edm::Event &iEvent, edm::EventSetup const &setup) const {
  using namespace edm;

  //get vertices
  Handle<reco::VertexCollection> vertexCollection;
  std::unique_ptr<Run3ScoutingVertexCollection> outVertices(new Run3ScoutingVertexCollection());
  if (iEvent.getByToken(vertexCollection_, vertexCollection)) {
    for (auto &vtx : *vertexCollection) {
        outVertices->emplace_back(MiniFloatConverter::reduceMantissaToNbitsRounding(vtx.x(), mantissaPrecision),
                                  MiniFloatConverter::reduceMantissaToNbitsRounding(vtx.y(), mantissaPrecision),
                                  MiniFloatConverter::reduceMantissaToNbitsRounding(vtx.z(), mantissaPrecision),
                                  MiniFloatConverter::reduceMantissaToNbitsRounding(vtx.zError(), mantissaPrecision),
                                  MiniFloatConverter::reduceMantissaToNbitsRounding(vtx.xError(), mantissaPrecision),
                                  MiniFloatConverter::reduceMantissaToNbitsRounding(vtx.yError(), mantissaPrecision),
                                  vtx.tracksSize(),
                                  MiniFloatConverter::reduceMantissaToNbitsRounding(vtx.chi2(), mantissaPrecision),
                                  vtx.ndof(),
                                  vtx.isValid());
    }
  }

  //get rho
  Handle<double> rho;
  std::unique_ptr<double> outRho(new double(-999));
  if (iEvent.getByToken(rho_, rho)) {
    outRho = std::make_unique<double>(*rho);
  }

  //get MET
  Handle<reco::PFMETCollection> metCollection;
  std::unique_ptr<double> outMetPt(new double(-999));
  std::unique_ptr<double> outMetPhi(new double(-999));
  if (doMet && iEvent.getByToken(metCollection_, metCollection)) {
    outMetPt = std::make_unique<double>(metCollection->front().pt());
    outMetPhi = std::make_unique<double>(metCollection->front().phi());
  }

  // get PF candidates with ParticleNet features
  edm::ESHandle<TransientTrackBuilder> trackBuilder = setup.getHandle(trackBuilderToken);
  TrackInfoBuilder trkinfo(trackBuilder);

  auto useTrackProperties = [&](const reco::PFCandidate *reco_cand) {
      const auto *trk = reco_cand->bestTrack();
      return trk != nullptr;
  };

  //get PF candidates
  Handle<reco::PFCandidateCollection> pfCandidateCollection;
  std::unique_ptr<Run3ScoutingParticleV2Collection> outPFCandidates(new Run3ScoutingParticleV2Collection());
  if (doCandidates && iEvent.getByToken(pfCandidateCollection_, pfCandidateCollection)) {
    for (auto &cand : *pfCandidateCollection) {
      if (cand.pt() > pfCandidatePtCut && std::abs(cand.eta()) < pfCandidateEtaCut) {

        /*

        int vertex_index = -1;
        int index_counter = 0;
        double dr2 = 0.0001;
        for (auto &vtx : *outVertices) {
          double tmp_dr2 = pow(vtx.x() - cand.vx(), 2) + pow(vtx.y() - cand.vy(), 2) + pow(vtx.z() - cand.vz(), 2);
          if (tmp_dr2 < dr2) {
            dr2 = tmp_dr2;
            vertex_index = index_counter;
          }
          if (dr2 == 0.0)
            break;
          ++index_counter;
        }

        outPFCandidates->emplace_back(MiniFloatConverter::reduceMantissaToNbitsRounding(cand.pt(), mantissaPrecision),
                                      MiniFloatConverter::reduceMantissaToNbitsRounding(cand.eta(), mantissaPrecision),
                                      MiniFloatConverter::reduceMantissaToNbitsRounding(cand.phi(), mantissaPrecision),
                                      MiniFloatConverter::reduceMantissaToNbitsRounding(cand.mass(), mantissaPrecision),
                                      cand.pdgId(),
                                      vertex_index);
        */

        const auto *reco_cand = dynamic_cast<const reco::PFCandidate*> (&cand);
        const reco::Vertex *pv = &vertexCollection->at(0);
        const auto *trk = cand.bestTrack();
        float dz = trk ? trk->dz(pv->position()) : 0;
        float dxy = trk ? trk->dxy(pv->position()) : 0;

        outPFCandidates->emplace_back(
            MiniFloatConverter::reduceMantissaToNbitsRounding(cand.pdgId(), mantissaPrecision),
            MiniFloatConverter::reduceMantissaToNbitsRounding(cand.mass(), mantissaPrecision),
            MiniFloatConverter::reduceMantissaToNbitsRounding(cand.pt(), mantissaPrecision),
            MiniFloatConverter::reduceMantissaToNbitsRounding(cand.eta(), mantissaPrecision),
            MiniFloatConverter::reduceMantissaToNbitsRounding(cand.phi(), mantissaPrecision),
            trk ? MiniFloatConverter::reduceMantissaToNbitsRounding(std::floor(trk->normalizedChi2()), mantissaPrecision) : MiniFloatConverter::reduceMantissaToNbitsRounding(999, mantissaPrecision), // normchi2,
            MiniFloatConverter::reduceMantissaToNbitsRounding(dz, mantissaPrecision),
            MiniFloatConverter::reduceMantissaToNbitsRounding(dxy, mantissaPrecision),
            trk ? MiniFloatConverter::reduceMantissaToNbitsRounding(dz / trk->dzError(), mantissaPrecision) : 0, // dzsig
            trk ? MiniFloatConverter::reduceMantissaToNbitsRounding(dxy / trk->dxyError(), mantissaPrecision) : 0, // dxysig
            useTrackProperties(reco_cand) ? lost_inner_hits_from_pfcand(*reco_cand) : 0, // lostInnerHits
            useTrackProperties(reco_cand) ? quality_from_pfcand(*reco_cand) : 0, // quality
            cand.charge(),
            trk ? MiniFloatConverter::reduceMantissaToNbitsRounding(trk->pt(), mantissaPrecision) : 0,
            trk ? MiniFloatConverter::reduceMantissaToNbitsRounding(trk->eta(), mantissaPrecision) : 0,
            trk ? MiniFloatConverter::reduceMantissaToNbitsRounding(trk->phi(), mantissaPrecision) : 0
            );

      }
    }
  }

  //get PF jets
  Handle<reco::PFJetCollection> pfJetCollection;
  std::unique_ptr<Run3ScoutingPFJetCollection> outPFJets(new Run3ScoutingPFJetCollection());
  if (iEvent.getByToken(pfJetCollection_, pfJetCollection)) {
    //get PF jet tags
    Handle<reco::JetTagCollection> pfJetTagCollection;
    bool haveJetTags = false;
    if (doJetTags && iEvent.getByToken(pfJetTagCollection_, pfJetTagCollection)) {
      haveJetTags = true;
    }

    for (auto &jet : *pfJetCollection) {
      if (jet.pt() < pfJetPtCut || std::abs(jet.eta()) > pfJetEtaCut)
        continue;
      //find the jet tag corresponding to the jet
      float tagValue = -20;
      float minDR2 = 0.01;
      if (haveJetTags) {
        for (auto &tag : *pfJetTagCollection) {
          float dR2 = reco::deltaR2(jet, *(tag.first));
          if (dR2 < minDR2) {
            minDR2 = dR2;
            tagValue = tag.second;
          }
        }
      }
      //get the PF constituents of the jet
      std::vector<int> candIndices;
      /*
      if (doCandidates) {
        for (auto &cand : jet.getPFConstituents()) {
          if (cand->pt() > pfCandidatePtCut && std::abs(cand->eta()) < pfCandidateEtaCut) {
            //search for the candidate in the collection
            float minDR2 = 0.0001;
            int matchIndex = -1;
            int outIndex = 0;
            for (auto &outCand : *outPFCandidates) {
              float dR2 = pow(cand->eta() - outCand.eta(), 2) + pow(cand->phi() - outCand.phi(), 2);
              if (dR2 < minDR2) {
                minDR2 = dR2;
                matchIndex = outIndex;
              }
              if (minDR2 == 0) {
                break;
              }
              outIndex++;
            }
            candIndices.push_back(matchIndex);
          }
        }
      }
      */
      outPFJets->emplace_back(MiniFloatConverter::reduceMantissaToNbitsRounding(jet.pt(), mantissaPrecision),
                              MiniFloatConverter::reduceMantissaToNbitsRounding(jet.eta(), mantissaPrecision),
                              MiniFloatConverter::reduceMantissaToNbitsRounding(jet.phi(), mantissaPrecision),
                              MiniFloatConverter::reduceMantissaToNbitsRounding(jet.mass(), mantissaPrecision),
                              MiniFloatConverter::reduceMantissaToNbitsRounding(jet.jetArea(), mantissaPrecision),
                              MiniFloatConverter::reduceMantissaToNbitsRounding(jet.chargedHadronEnergy(), mantissaPrecision),
                              MiniFloatConverter::reduceMantissaToNbitsRounding(jet.neutralHadronEnergy(), mantissaPrecision),
                              MiniFloatConverter::reduceMantissaToNbitsRounding(jet.photonEnergy(), mantissaPrecision),
                              MiniFloatConverter::reduceMantissaToNbitsRounding(jet.electronEnergy(), mantissaPrecision),
                              MiniFloatConverter::reduceMantissaToNbitsRounding(jet.muonEnergy(), mantissaPrecision),
                              MiniFloatConverter::reduceMantissaToNbitsRounding(jet.HFHadronEnergy(), mantissaPrecision),
                              MiniFloatConverter::reduceMantissaToNbitsRounding(jet.HFEMEnergy(), mantissaPrecision),
                              jet.chargedHadronMultiplicity(),
                              jet.neutralHadronMultiplicity(),
                              jet.photonMultiplicity(),
                              jet.electronMultiplicity(),
                              jet.muonMultiplicity(),
                              jet.HFHadronMultiplicity(),
                              jet.HFEMMultiplicity(),
                              MiniFloatConverter::reduceMantissaToNbitsRounding(jet.hoEnergy(), mantissaPrecision),
                              tagValue,
                              0.0,
                              candIndices);
    }
  }

  //put output
  iEvent.put(std::move(outPFCandidates));
  iEvent.put(std::move(outPFJets));
  iEvent.put(std::move(outRho), "rho");
  iEvent.put(std::move(outMetPt), "pfMetPt");
  iEvent.put(std::move(outMetPhi), "pfMetPhi");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HLTScoutingPFProducerV2::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pfJetCollection", edm::InputTag("hltAK4PFJets"));
  desc.add<edm::InputTag>("pfJetTagCollection", edm::InputTag("hltDeepCombinedSecondaryVertexBJetTagsPF"));
  desc.add<edm::InputTag>("pfCandidateCollection", edm::InputTag("hltParticleFlow"));
  desc.add<edm::InputTag>("vertexCollection", edm::InputTag("hltPixelVertices"));
  desc.add<edm::InputTag>("metCollection", edm::InputTag("hltPFMETProducer"));
  desc.add<edm::InputTag>("rho", edm::InputTag("hltFixedGridRhoFastjetAll"));
  desc.add<double>("pfJetPtCut", 20.0);
  desc.add<double>("pfJetEtaCut", 3.0);
  desc.add<double>("pfCandidatePtCut", 0.6);
  desc.add<double>("pfCandidateEtaCut", 5.0);
  desc.add<int>("mantissaPrecision", 10)->setComment("default float16, change to 23 for float32");
  desc.add<bool>("doJetTags", true);
  desc.add<bool>("doCandidates", true);
  desc.add<bool>("doMet", true);
  descriptions.add("hltScoutingPFProducerV2", desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(HLTScoutingPFProducerV2);
