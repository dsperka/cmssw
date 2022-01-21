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
//#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Math/interface/libminifloat.h"

class HLTScoutingPackedCandidateProducer : public edm::global::EDProducer<> {
public:
  explicit HLTScoutingPackedCandidateProducer(const edm::ParameterSet &);
  ~HLTScoutingPackedCandidateProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    // sorting of cands to maximize the zlib compression
    static bool candsOrdering(pat::PackedCandidate const &i, pat::PackedCandidate const &j) {
        if (std::abs(i.charge()) == std::abs(j.charge())) {
            if (i.charge() != 0) {
                if (i.hasTrackDetails() and !j.hasTrackDetails())
                    return true;
                if (!i.hasTrackDetails() and j.hasTrackDetails())
                    return false;
                if (i.covarianceSchema() > j.covarianceSchema())
                    return true;
                if (i.covarianceSchema() < j.covarianceSchema())
                    return false;
            }
            if (i.vertexRef() == j.vertexRef())
                return i.eta() > j.eta();
            else
                return i.vertexRef().key() < j.vertexRef().key();
        }
        return std::abs(i.charge()) > std::abs(j.charge());
    }

    template <typename T>
    static std::vector<size_t> sort_indexes(const std::vector<T> &v) {
        std::vector<size_t> idx(v.size());
        for (size_t i = 0; i != idx.size(); ++i)
            idx[i] = i;
        std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return candsOrdering(v[i1], v[i2]); });
        return idx;
    }

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


    const double minPtForChargedHadronProperties_;
    const double minPtForTrackProperties_;
    const double minPtForLowQualityTrackProperties_;
    const int covarianceVersion_;
    const std::vector<int> covariancePackingSchemas_;

    const std::vector<int> pfCandidateTypesForHcalDepth_;
    const bool storeHcalDepthEndcapOnly_;

};

//
// constructors and destructor
//
HLTScoutingPackedCandidateProducer::HLTScoutingPackedCandidateProducer(const edm::ParameterSet &iConfig)
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
    minPtForChargedHadronProperties_(iConfig.getParameter<double>("minPtForChargedHadronProperties")),
    minPtForTrackProperties_(iConfig.getParameter<double>("minPtForTrackProperties")),
    minPtForLowQualityTrackProperties_(iConfig.getParameter<double>("minPtForLowQualityTrackProperties")),
    covarianceVersion_(iConfig.getParameter<int>("covarianceVersion")),
    covariancePackingSchemas_(iConfig.getParameter<std::vector<int>>("covariancePackingSchemas")),
    pfCandidateTypesForHcalDepth_(iConfig.getParameter<std::vector<int>>("pfCandidateTypesForHcalDepth")),
    storeHcalDepthEndcapOnly_(iConfig.getParameter<bool>("storeHcalDepthEndcapOnly"))
 {
  //register products
  produces<Run3ScoutingPFJetCollection>();
  //produces<Run3ScoutingParticleCollection>();
  produces<std::vector<pat::PackedCandidate>>();
  produces<double>("rho");
  produces<double>("pfMetPt");
  produces<double>("pfMetPhi");
}

HLTScoutingPackedCandidateProducer::~HLTScoutingPackedCandidateProducer() = default;

// ------------ method called to produce the data  ------------
void HLTScoutingPackedCandidateProducer::produce(edm::StreamID sid, edm::Event &iEvent, edm::EventSetup const &setup) const {
  using namespace edm;

  //get vertices
  Handle<reco::VertexCollection> vertexCollection;
  std::unique_ptr<Run3ScoutingVertexCollection> outVertices(new Run3ScoutingVertexCollection());

  if (iEvent.getByToken(vertexCollection_, vertexCollection)) {
    for (auto &vtx : *vertexCollection) {
      outVertices->emplace_back(vtx.x(),
                                vtx.y(),
                                vtx.z(),
                                vtx.zError(),
                                vtx.xError(),
                                vtx.yError(),
                                vtx.tracksSize(),
                                vtx.chi2(),
                                vtx.ndof(),
                                vtx.isValid());
    }
  }

  reco::VertexRef PV(vertexCollection.id());
  reco::VertexRefProd PVRefProd(vertexCollection);
  math::XYZPoint PVpos;


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

  //get PF candidates
  Handle<reco::PFCandidateCollection> pfCandidateCollection;
  //std::unique_ptr<Run3ScoutingParticleCollection> outPFCandidates(new Run3ScoutingParticleCollection());
  auto outPFCandidates = std::make_unique<std::vector<pat::PackedCandidate>>();
  auto outPFCandidatesSorted = std::make_unique<std::vector<pat::PackedCandidate>>();

  if (doCandidates && iEvent.getByToken(pfCandidateCollection_, pfCandidateCollection)) {


    //for (auto &cand : *pfCandidateCollection) {        
    for (unsigned int ic = 0, nc = pfCandidateCollection->size(); ic < nc; ++ic) {
      const reco::PFCandidate &cand = (*pfCandidateCollection)[ic];
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


        const reco::Track *ctrack = nullptr;
        if ((abs(cand.pdgId()) == 11 || cand.pdgId() == 22) && cand.gsfTrackRef().isNonnull()) {
            ctrack = &*cand.gsfTrackRef();
        } else if (cand.trackRef().isNonnull()) {
            ctrack = &*cand.trackRef();
        }
        if (ctrack) {
            float dist = 1e99;
            int pvi = -1;

            for (size_t ii = 0; ii < vertexCollection->size(); ii++) {
                float dz = std::abs(ctrack->dz(((*vertexCollection)[ii]).position()));
                if (dz < dist) {
                    pvi = ii;
                    dist = dz;
                }
            }

            PV = reco::VertexRef(vertexCollection, pvi);
            math::XYZPoint vtx = cand.vertex();
            pat::PackedCandidate::LostInnerHits lostHits = pat::PackedCandidate::noLostInnerHits;
            
            /*
            const reco::VertexRef &PVOrig = associatedPV[reco::CandidatePtr(pfCandidateCollection, ic)];
            if (PVOrig.isNonnull())
                PV = reco::VertexRef(vertexCollection,
                                     PVOrig.key());  // WARNING: assume the PV slimmer is keeping same order

            int quality = associationQuality[reco::CandidatePtr(pfCandidateCollection, ic)];
            */

            vtx = ctrack->referencePoint();
            float ptTrk = ctrack->pt();
            float etaAtVtx = ctrack->eta();
            float phiAtVtx = ctrack->phi();

            int nlost = ctrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
            if (nlost == 0) {
                if (ctrack->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1)) {
                    lostHits = pat::PackedCandidate::validHitInFirstPixelBarrelLayer;
                }
            } else {
                lostHits = (nlost == 1 ? pat::PackedCandidate::oneLostInnerHit : pat::PackedCandidate::moreLostInnerHits);
            }

            outPFCandidates->push_back(
                pat::PackedCandidate(cand.polarP4(), vtx, ptTrk, etaAtVtx, phiAtVtx, cand.pdgId(), PVRefProd, PV.key()));
            
            
            outPFCandidates->back().setCovarianceVersion(covarianceVersion_);
            
            /*
            const static int qualityMap[8] = {1, 0, 1, 1, 4, 4, 5, 6};
            outPFCandidates->back().setAssociationQuality(pat::PackedCandidate::PVAssociationQuality(qualityMap[quality]));
            if (cand.trackRef().isNonnull() && PVOrig.isNonnull() && PVOrig->trackWeight(cand.trackRef()) > 0.5 &&
                quality == 7) {
                outPFCandidates->back().setAssociationQuality(pat::PackedCandidate::UsedInFitTight);
            }
            */

            // properties of the best track
            outPFCandidates->back().setLostInnerHits(lostHits);
            if (outPFCandidates->back().pt() > minPtForTrackProperties_ || outPFCandidates->back().ptTrk() > minPtForTrackProperties_) {
                outPFCandidates->back().setFirstHit(ctrack->hitPattern().getHitPattern(reco::HitPattern::TRACK_HITS, 0));
                if (abs(outPFCandidates->back().pdgId()) == 22) {
                    outPFCandidates->back().setTrackProperties(*ctrack, covariancePackingSchemas_[4], covarianceVersion_);
                } else {
                    if (ctrack->hitPattern().numberOfValidPixelHits() > 0) {
                        outPFCandidates->back().setTrackProperties(*ctrack,
                                                           covariancePackingSchemas_[0],
                                                           covarianceVersion_);  // high quality
                    } else {
                        outPFCandidates->back().setTrackProperties(*ctrack, covariancePackingSchemas_[1], covarianceVersion_);
                    }
                }
                // outPFCandidates->back().setTrackProperties(*ctrack,tsos.curvilinearError());
            } else {
                if (outPFCandidates->back().pt() > minPtForLowQualityTrackProperties_) {
                    if (ctrack->hitPattern().numberOfValidPixelHits() > 0)
                        outPFCandidates->back().setTrackProperties(*ctrack,
                                                           covariancePackingSchemas_[2],
                                                           covarianceVersion_);  // low quality, with pixels
                    else
                        outPFCandidates->back().setTrackProperties(*ctrack,
                                                           covariancePackingSchemas_[3],
                                                           covarianceVersion_);  // low quality, without pixels
                }
            }

            // these things are always for the CKF track
            outPFCandidates->back().setTrackHighPurity(cand.trackRef().isNonnull() &&
                                               cand.trackRef()->quality(reco::Track::highPurity));
            if (cand.muonRef().isNonnull()) {
                outPFCandidates->back().setMuonID(cand.muonRef()->isStandAloneMuon(), cand.muonRef()->isGlobalMuon());
            }
        } else {
            if (!vertexCollection->empty()) {
                PV = reco::VertexRef(vertexCollection, 0);
                PVpos = PV->position();
            }

            outPFCandidates->push_back(pat::PackedCandidate(
                                   cand.polarP4(), PVpos, cand.pt(), cand.eta(), cand.phi(), cand.pdgId(), PVRefProd, PV.key()));
            outPFCandidates->back().setAssociationQuality(
                pat::PackedCandidate::PVAssociationQuality(pat::PackedCandidate::UsedInFitTight));
        }

        // neutrals and isolated charged hadrons

        bool isIsolatedChargedHadron = false;
        
        /*
        if (storeChargedHadronIsolation_) {
            const edm::ValueMap<bool> &chargedHadronIsolation = *(chargedHadronIsolationHandle.product());
      isIsolatedChargedHadron =
          ((cand.pt() > minPtForChargedHadronProperties_) && (chargedHadronIsolation[reco::PFCandidateRef(cands, ic)]));
      outPFCandidates->back().setIsIsolatedChargedHadron(isIsolatedChargedHadron);
        }
        */

        if (abs(cand.pdgId()) == 1 || abs(cand.pdgId()) == 130) {
            outPFCandidates->back().setHcalFraction(cand.hcalEnergy() / (cand.ecalEnergy() + cand.hcalEnergy()));
        } else if ((cand.charge() || abs(cand.pdgId()) == 22) && cand.pt() > 0.5) {
            outPFCandidates->back().setHcalFraction(cand.hcalEnergy() / (cand.ecalEnergy() + cand.hcalEnergy()));
            outPFCandidates->back().setCaloFraction((cand.hcalEnergy() + cand.ecalEnergy()) / cand.energy());
        } else {
            outPFCandidates->back().setHcalFraction(0);
            outPFCandidates->back().setCaloFraction(0);
        }

        if (isIsolatedChargedHadron) {
            outPFCandidates->back().setRawCaloFraction((cand.rawEcalEnergy() + cand.rawHcalEnergy()) / cand.energy());
            outPFCandidates->back().setRawHcalFraction(cand.rawHcalEnergy() / (cand.rawEcalEnergy() + cand.rawHcalEnergy()));
        } else {
            outPFCandidates->back().setRawCaloFraction(0);
            outPFCandidates->back().setRawHcalFraction(0);
        }

        /*
        std::vector<float> dummyVector;
        dummyVector.clear();
        pat::HcalDepthEnergyFractions hcalDepthEFrac(dummyVector);

        // storing HcalDepthEnergyFraction information
        if (std::find(pfCandidateTypesForHcalDepth_.begin(), pfCandidateTypesForHcalDepth_.end(), abs(cand.pdgId())) !=
            pfCandidateTypesForHcalDepth_.end()) {
            if (!storeHcalDepthEndcapOnly_ ||
                fabs(outPFCandidates->back().eta()) > 1.3) {  // storeHcalDepthEndcapOnly_==false -> store all eta of
                // selected PF types, if true, only |eta|>1.3 of selected
                // PF types will be stored
                std::vector<float> hcalDepthEnergyFractionTmp(cand.hcalDepthEnergyFractions().begin(),
                                                              cand.hcalDepthEnergyFractions().end());
                hcalDepthEFrac.reset(hcalDepthEnergyFractionTmp);
            }
        }
        hcalDepthEnergyFractions.push_back(hcalDepthEFrac);
        */

        // specifically this is the PFLinker requirements to apply the e/gamma
        // regression
        if (cand.particleId() == reco::PFCandidate::e ||
            (cand.particleId() == reco::PFCandidate::gamma && cand.mva_nothing_gamma() > 0.)) {
            outPFCandidates->back().setGoodEgamma();
        }

        /*
        mapping[ic] = ic;  // trivial at the moment!
        if (cand.trackRef().isNonnull() && cand.trackRef().id() == TKOrigs.id()) {
            mappingTk[cand.trackRef().key()] = ic;
        }
        */
      }

    }

    std::vector<size_t> order = sort_indexes(*outPFCandidates);
    std::vector<size_t> reverseOrder(order.size());
    for (size_t i = 0, nc = outPFCandidates->size(); i < nc; i++) {
        outPFCandidatesSorted->push_back((*outPFCandidates)[order[i]]);
        reverseOrder[order[i]] = i;
        //mappingReverse[order[i]] = i;
        //hcalDepthEnergyFractions_Ordered.push_back(hcalDepthEnergyFractions[order[i]]);
    }
    
    // Fix track association for sorted candidates
    /*
      for (size_t i = 0, ntk = mappingTk.size(); i < ntk; i++) {
      if (mappingTk[i] >= 0)
      mappingTk[i] = reverseOrder[mappingTk[i]];
      }
    */
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
      outPFJets->emplace_back(jet.pt(),
                              jet.eta(),
                              jet.phi(),
                              jet.mass(),
                              jet.jetArea(),
                              jet.chargedHadronEnergy(),
                              jet.neutralHadronEnergy(),
                              jet.photonEnergy(),
                              jet.electronEnergy(),
                              jet.muonEnergy(),
                              jet.HFHadronEnergy(),
                              jet.HFEMEnergy(),
                              jet.chargedHadronMultiplicity(),
                              jet.neutralHadronMultiplicity(),
                              jet.photonMultiplicity(),
                              jet.electronMultiplicity(),
                              jet.muonMultiplicity(),
                              jet.HFHadronMultiplicity(),
                              jet.HFEMMultiplicity(),
                              jet.hoEnergy(),
                              tagValue,
                              0.0,
                              candIndices);
    }
  }

  //put output
  //iEvent.put(std::move(outPFCandidates));
  //edm::OrphanHandle<pat::PackedCandidateCollection> oh = iEvent.put(std::move(outPFCandidatesSorted));
  iEvent.put(std::move(outPFCandidatesSorted));
  iEvent.put(std::move(outPFJets));
  iEvent.put(std::move(outRho), "rho");
  iEvent.put(std::move(outMetPt), "pfMetPt");
  iEvent.put(std::move(outMetPhi), "pfMetPhi");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HLTScoutingPackedCandidateProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
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

  desc.add<double>("minPtForChargedHadronProperties",3.0);
  desc.add<double>("minPtForTrackProperties",0.3);
  desc.add<double>("minPtForLowQualityTrackProperties",0.3);
  desc.add<int>("covarianceVersion",0);

  std::vector<int> covariancePackingSchemas = {8,264,520,776,0};
  desc.add<std::vector<int> >("covariancePackingSchemas",covariancePackingSchemas);

  std::vector<int> pfCandidateTypesForHcalDepth = {};
  desc.add<std::vector<int> >("pfCandidateTypesForHcalDepth",pfCandidateTypesForHcalDepth);
  desc.add<bool>("storeHcalDepthEndcapOnly",false);

  descriptions.add("hltScoutingPackedCandidateProducer", desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(HLTScoutingPackedCandidateProducer);
