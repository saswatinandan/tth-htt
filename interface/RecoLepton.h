#ifndef tthAnalysis_HiggsToTauTau_RecoLepton_h
#define tthAnalysis_HiggsToTauTau_RecoLepton_h

#include "tthAnalysis/HiggsToTauTau/interface/ChargedParticle.h" // ChargedParticle

#include <memory> // std::shared_ptr<>
#include <map> // std::map<,>

// forward declarations
class GenLepton;
class GenHadTau;
class GenPhoton;
class GenJet;

enum class Btag;

class RecoLepton
  : public ChargedParticle
{
public:
  RecoLepton() = default;
  RecoLepton(const ChargedParticle & lepton,
             Double_t dxy,
             Double_t dz,
             Double_t relIso,
             Double_t pfRelIso04All,
             Double_t miniRelIsoCharged,
             Double_t miniRelIsoNeutral,
             Double_t sip3d,
             Double_t mvaRawTTH,
             Double_t jetPtRatio,
             Double_t jetPtRel,
             Int_t    jetNDauChargedMVASel,
             Int_t    tightCharge,
             UInt_t   filterBits,
             Int_t    jetIdx,
             UChar_t  genPartFlav,
             Int_t    genMatchIdx);

  virtual ~RecoLepton();

  /**
   * @brief Set flags indicating whether or not lepton passes loose, fakeable and/or tight selection criteria
   */
  void set_isCMSPOG() const;
  void set_isLoose() const;
  void set_isFakeable() const;
  void set_isTight() const;

  /**
   * @brief Set links to generator level particles (matched by dR)
   */
  void set_genLepton(const GenLepton * genLepton);
  void set_genHadTau(const GenHadTau * genHadTau);
  void set_genPhoton(const GenPhoton * genPhoton);
  void set_genJet(const GenJet * genJet);

  /**
   * @brief Checks whether a given lepton is an electron by its PDG id
   * @return True if it is an electron; false otherwise
   */
  virtual bool
  is_electron() const;

  /**
   * @brief Checks whether a given lepton is a muon by its PDG id
   * @return True if it is a muon; false otherwise
   */
  virtual bool
  is_muon() const;

  /**
   * @brief Funtions to access pT and four-momenta
   * @return Values of data-members
   *
   * CV: note that leptons have different pT values, 
   *     depending on whether or not they pass the "tight" lepton selection criteria
   *    (cf. lines 289-291 of AN-2017/029 v5)
   *     The factor 0.85 has been updated to 0.90 for the analysis of 2017 data.
   *     The different pT definitions are as follows:
   *       - The lepton_pt() and pt() functions always return the pT of the reco lepton
   *       - The cone_pt() function returns the pT of the reco lepton if the lepton passes the "tight" lepton MVA
   *         (and if the lepton is a muon, it is required to pass Medium muon ID definition), otherwise 0.90*jet pT is returned
   *       - The assocJet_pt() function always returns 0.90*jet pT
   *     Giovanni clarified via email on 02/23/2018 that the cone_pt() function is to be used everywhere in the analysis:
   *       - for sorting the leptons by decreasing pT
   *       - for evaluating the event level BDTs
   *       - for applying fake-rates in the analysis 
   *        (Note: fake-rates are only applied to leptons passing "fakeable" and failing "tight" lepton selection,
   *               for leptons passing the "tight" selection criteria the fake-rate weight is unity)
   *       - for applying pT cuts in the analysis (e.g. 25/15 GeV for leading/subleading lepton)
   *     The assocJet_pt() function is used only for the purpose of filling the numerator and denominator histograms
   *     in the code that computes the fake-rates (using the assocJet_pt() instead of the cone_pt() function guarantees
   *     that the same pT definition is used in the numerator and denominator histograms, 
   *     i.e. regardless of whether leptons pass or fail the "tight" lepton selection criteria.
   */

  Double_t
  lepton_pt() const;

  const Particle::LorentzVector &
  lepton_p4() const;

  virtual Double_t
  cone_pt() const;

  virtual const Particle::LorentzVector &
  cone_p4() const;

  Double_t
  assocJet_pt() const;

  const Particle::LorentzVector &
  assocJet_p4() const;

  /**
   * @brief Funtions to access other data-members
   * @return Values of data-members
   */

  Double_t dxy() const;
  Double_t dz() const;
  Double_t relIso() const;
  Double_t pfRelIso04All() const;
  Double_t miniRelIsoCharged() const;
  Double_t miniRelIsoNeutral() const;
  Double_t sip3d() const;
  Double_t mvaRawTTH() const;
  Double_t mvaRawTTH_cut() const;
  Double_t jetPtRatio() const;
  Double_t jetRelIso() const;
  Double_t jetPtRel() const;
  Double_t jetBtagCSV(bool doAssoc = false) const;
  Double_t jetBtagCSV(Btag btag, bool doAssoc = false) const;
  Int_t jetNDauChargedMVASel() const;
  Int_t tightCharge() const;
  UInt_t filterBits() const;
  Int_t jetIdx() const;
  UChar_t genPartFlav() const;
  Int_t genMatchIdx() const;

  const GenLepton * genLepton() const;
  const GenHadTau * genHadTau() const;
  const GenPhoton * genPhoton() const;
  const GenJet * genJet() const;

  bool hasJetBtagCSV(Btag btag, bool doAssoc = false) const;

  bool isGenMatched(bool requireChargeMatch) const;
  bool hasAnyGenMatch() const;

  bool isCMSPOG() const;
  bool isLoose() const;
  bool isFakeable() const;
  bool isTight() const;

  void set_mvaRawTTH_cut(Double_t mvaRawTTH_cut);

  virtual void
  set_p4(const Particle::LorentzVector & p4) override;

  virtual void
  set_ptEtaPhiMass(Double_t pt,
                   Double_t eta,
                   Double_t phi,
                   Double_t mass) override;

  friend class RecoMuonReader;
  friend class RecoElectronReader;

protected:
//--- common observables for electrons and muons
  Double_t dxy_;                ///< d_{xy}, distance in the transverse plane w.r.t PV
  Double_t dz_;                 ///< d_{z}, distance on the z axis w.r.t PV
  Double_t relIso_;             ///< relative mini-isolation
  Double_t pfRelIso04All_;      ///< PF relative isolation dR=0.3, charged component
  Double_t miniRelIsoCharged_;  ///< relative charged mini-isolation
  Double_t miniRelIsoNeutral_;  ///< relative neutral mini-isolation (PU corrected)
  Double_t sip3d_;              ///< significance of IP
  Double_t mvaRawTTH_;          ///< raw output of lepton MVA of ttH multilepton analysis
  Double_t jetPtRatio_;         ///< ratio of lepton pT to pT of nearby jet
  Double_t jetPtRel_;           ///< perpendicular component of the distance vector between lepton and its jet pT vectors
  Int_t jetNDauChargedMVASel_;  ///< number of charged constituents in the nearest jet
  Int_t tightCharge_;           ///< Flag indicating if lepton passes (>= 2) or fails (< 2) tight charge requirement
  UInt_t filterBits_;           ///< bitmask of matching with trigger objects
  Int_t jetIdx_;                ///< index of jet from initial jet collection that the lepton is constituent of (-1 if no match)
  UChar_t genPartFlav_;         ///< generator-level parton flavor
  Int_t genMatchIdx_;           ///< index to matched gen particle (-1 if no match)
  Double_t mvaRawTTH_cut_;      ///< cut on prompt lepton MVA score

  std::map<Btag, Double_t> jetBtagCSVs_; ///< CSV b-tagging discriminator values of nearby jet as used in prompt lepton MVA
  std::map<Btag, Double_t> assocJetBtagCSVs_; ///< CSV b-tagging discriminator values of nearby jet found via jetIdx branch

  Double_t assocJet_pt_;
  Particle::LorentzVector assocJet_p4_;

//--- matching to generator level particles
  std::shared_ptr<const GenLepton> genLepton_;
  std::shared_ptr<const GenHadTau> genHadTau_;
  std::shared_ptr<const GenPhoton> genPhoton_;
  std::shared_ptr<const GenJet> genJet_;

//--- flags indicating whether or not lepton passes CMS POG ID, loose, fakeable and/or tight selection criteria
  mutable bool isCMSPOG_;
  mutable bool isLoose_;
  mutable bool isFakeable_;
  mutable bool isTight_;

  void
  set_assocJet_p4();

  static Double_t
  get_assocJet_pt(Double_t reco_pt,
                  Double_t jetPtRatio);
};

std::ostream &
operator<<(std::ostream & stream,
           const RecoLepton & lepton);

#endif // tthAnalysis_HiggsToTauTau_RecoLepton_h
