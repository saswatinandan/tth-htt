#ifndef tthAnalysis_HiggsToTauTau_RecoElectronCollectionSelectorTight_h
#define tthAnalysis_HiggsToTauTau_RecoElectronCollectionSelectorTight_h

#include "tthAnalysis/HiggsToTauTau/interface/ParticleCollectionSelector.h" // ParticleCollectionSelector
#include "tthAnalysis/HiggsToTauTau/interface/RecoElectron.h"               // RecoElectron
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h"       // Era

class RecoElectronSelectorTight
{
public:
  RecoElectronSelectorTight(Era era,
                            int index = -1,
                            bool debug = false,
                            bool set_selection_flags = true);

  ~RecoElectronSelectorTight() {}

  /**
   * @brief Set cut thresholds
   */
  // enable/disable cuts on electron ID variables to mimic electron ID cuts applied by single electron trigger 
  void enable_offline_e_trigger_cuts();
  void disable_offline_e_trigger_cuts();

  // enable/disable photon conversion veto
  void enable_conversionVeto();
  void disable_conversionVeto();

  void invert_max_nLostHits(Int_t min_nLostHits_fornLostHitsInversion); // for conversion bkg CR
  void invert_conversionVeto();

  void set_min_lepton_pt(double min_lepton_pt);
  void set_min_cone_pt(double min_cone_pt);
  void set_max_absEta(double max_absEta);

  void set_selection_flags(bool selection_flag);

  void set_assocJetBtag(bool flag);

  /**
   * @brief Get cut thresholds
   */
  double get_min_lepton_pt() const;
  double get_min_cone_pt() const;
  double get_max_absEta() const;

  /**
   * @brief Check if electron given as function argument passes "tight" electron selection, defined in Table 13 of AN-2015/321
   * @return True if electron passes selection; false otherwise
   */
  bool operator()(const RecoElectron & electron) const;

  void print_selection_conditions();
  
protected:
  const Era era_;
  bool set_selection_flags_;
  bool apply_offline_e_trigger_cuts_;
  bool debug_;

  Double_t min_lepton_pt_;                   ///< lower cut threshold on reco::GSFElectron pT
  Double_t min_cone_pt_;                     ///< lower cut threshold on cone pT
  Double_t max_absEta_;                      ///< upper cut threshold on absolute value of eta
  const Double_t max_dxy_;                   ///< upper cut threshold on d_{xy}, distance in the transverse plane w.r.t PV
  const Double_t max_dz_;                    ///< upper cut threshold on d_{z}, distance on the z axis w.r.t PV
  const Double_t max_relIso_;                ///< upper cut threshold on relative isolation
  const Double_t max_sip3d_;                 ///< upper cut threshold on significance of IP
//--- define cuts that dependent on eta
//    format: central region (|eta| < 0.8) / transition region (0.8 < |eta| < 1.479) / forward region (|eta| > 1.479)
  const Double_t binning_absEta_;            ///< eta values separating transition and forward region
  const Double_t min_sigmaEtaEta_trig_;      ///< lower cut threshold on second shower moment in eta-direction
  const Double_t max_sigmaEtaEta_trig_;      ///< upper cut threshold on second shower moment in eta-direction
  const Double_t max_HoE_trig_;              ///< upper cut threshold on ratio of energy deposits in hadronic/electromagnetic section of calorimeter
  const Double_t min_OoEminusOoP_trig_;      ///< lower cut threshold on difference between calorimeter energy and track momentum (1/E - 1/P)
//-------------------------------------------------------------------------------
  const Double_t max_jetBtagCSV_;            ///< upper cut threshold on CSV b-tagging discriminator value of nearby jet
  const Int_t max_nLostHits_;                ///< upper cut threshold on lost hits in the innermost layer of the tracker (electrons with lost_hits equal to cut threshold pass)
  bool apply_conversionVeto_;                ///< apply (True) or do not apply (False) conversion veto
  bool invert_nLostHits_;                     /// for conversion bkg CR
  Int_t min_nLostHits_fornLostHitsInversion_; ///<fornLostHitsInversion: lower cut threshold on lost hits in the innermost layer of the tracker (electrons with lost_hits equal to cut threshold pass)
  bool invert_conversionVeto_;                ///< invert (True) or do not invert (False) conversion veto
//-------------------------------------------------------------------------------
  bool useAssocJetBtag_;                     ///< if true, use finalJets instead of updatedJets
};

class RecoElectronCollectionSelectorTight
  : public ParticleCollectionSelector<RecoElectron, RecoElectronSelectorTight>
{
public:
  explicit
  RecoElectronCollectionSelectorTight(Era era,
                                      int index = -1,
                                      bool debug = false,
                                      bool set_selection_flags = true);
  ~RecoElectronCollectionSelectorTight() {}

  // enable/disable cuts on electron ID variables to mimic electron ID cuts applied by single electron trigger 
  void enable_offline_e_trigger_cuts();
  void disable_offline_e_trigger_cuts();

  // enable/disable photon conversion veto
  void enable_conversionVeto();
  void disable_conversionVeto();

  void invert_max_nLostHits(Int_t min_nLostHits_fornLostHitsInversion); // for conversion bkg CR
  void invert_conversionVeto();

  void print_selection_conditions();
};

#endif // tthAnalysis_HiggsToTauTau_RecoElectronCollectionSelectorTight_h

