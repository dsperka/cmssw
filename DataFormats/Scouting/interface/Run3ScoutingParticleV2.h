#ifndef DataFormats_Run3ScoutingParticleV2_h
#define DataFormats_Run3ScoutingParticleV2_h

#include <vector>

//class for holding PF candidate information, for use in data scouting
//IMPORTANT: the content of this class should be changed only in backwards compatible ways!
class Run3ScoutingParticleV2 {
public:
    //constructor with values for all data fields
    Run3ScoutingParticleV2(int pdgId,
                                    float m,
                                    float pt,
                                    float eta,
                                    float phi,
                                    float normchi2,
                                    float dz,
                                    float dxy,
                                    float dzsig,
                                    float dxysig,
                                    float lostInnerHits,
                                    float quality,
                                    float charge,
                                    float trk_pt,
                                    float trk_eta,
                                    float trk_phi)
        : pdgId_(pdgId),
        m_(m),
        pt_(pt),
        eta_(eta),
        phi_(phi),
        normchi2_(normchi2),
        dz_(dz),
        dxy_(dxy),
        dzsig_(dzsig),
        dxysig_(dxysig),
        lostInnerHits_(lostInnerHits),
        quality_(quality),
        charge_(charge),
        trk_pt_(trk_pt),
        trk_eta_(trk_eta),
        trk_phi_(trk_phi) {}

    // default constractor
    Run3ScoutingParticleV2()
        : pdgId_(0),
        m_(0),
        pt_(0),
        eta_(0),
        phi_(0),
        normchi2_(0),
        dz_(0),
        dxy_(0),
        dzsig_(0),
        dxysig_(0),
        lostInnerHits_(0),
        quality_(0),
        charge_(0),
        trk_pt_(0),
        trk_eta_(0),
        trk_phi_(0) {}

    //accessor functions
    int pdgId() const { return pdgId_; }
    float m() const { return m_; }
    float pt() const { return pt_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    float normchi2() const { return normchi2_; }
    float dz() const { return dz_; }
    float dxy() const { return dxy_; }
    float dzsig() const { return dzsig_; }
    float dxysig() const { return dxysig_; }
    float lostInnerHits() const { return lostInnerHits_; }
    float quality() const { return quality_; }
    float charge() const { return charge_; }
    float trk_pt() const { return trk_pt_; }
    float trk_eta() const { return trk_eta_; }
    float trk_phi() const { return trk_phi_; }

private:
    int pdgId_;
    float m_;
    float pt_;
    float eta_;
    float phi_;
    float normchi2_;
    float dz_;
    float dxy_;
    float dzsig_;
    float dxysig_;
    float lostInnerHits_;
    float quality_;
    float charge_;
    float trk_pt_;
    float trk_eta_;
    float trk_phi_;
};

typedef std::vector<Run3ScoutingParticleV2> Run3ScoutingParticleV2Collection;

#endif
