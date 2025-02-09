#ifndef tthAnalysis_HiggsToTauTau_Data_to_MC_CorrectionInterface_Base_h
#define tthAnalysis_HiggsToTauTau_Data_to_MC_CorrectionInterface_Base_h

#include "tthAnalysis/HiggsToTauTau/interface/sysUncertOptions.h" // FRet, FRmt

#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // Era, pileupJetID

#include <FWCore/ParameterSet/interface/ParameterSet.h> // edm::ParameterSet

// forward declarations
class TFile;
class TauIDSFTool;
class lutWrapperBase;
enum class TauID;

class Data_to_MC_CorrectionInterface_Base
{
public:
  Data_to_MC_CorrectionInterface_Base(Era era, const edm::ParameterSet & cfg);
  virtual ~Data_to_MC_CorrectionInterface_Base();

  //-----------------------------------------------------------------------------
  // overwrite configuration parameters (needed by analyze_jetToTauFakeRate.cc)
  void
  setHadTauSelection(const std::string & hadTauSelection);
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // set leptons, taus, and jets
  // (to be called once per event, before calling any of the getSF.. functions)
  void
  setLeptons(const std::vector<const RecoLepton *> & leptons,
             bool requireChargeMatch = false);

  void
  setHadTaus(const std::vector<const RecoHadTau *> & hadTaus);

  void
  setJets(const std::vector<const RecoJet *> & jets);
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // data/MC correction for electron and muon trigger efficiency
  virtual double
  getSF_leptonTriggerEff(TriggerSFsys central_or_shift) const;
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // data/MC corrections for electron and muon identification and isolation efficiency,
  // including the cut on the ttH multilepton MVA
  double
  getSF_leptonID_and_Iso_loose(LeptonIDSFsys central_or_shift) const;

  double
  getSF_leptonID_and_Iso_looseToFakeable() const;

  double
  getSF_leptonID_and_Iso_tight_to_loose_woTightCharge(LeptonIDSFsys central_or_shift) const;

  double
  getSF_leptonID_and_Iso_tight_to_loose_wTightCharge(LeptonIDSFsys central_or_shift) const;
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // data/MC corrections for hadronic tau identification efficiency,
  // and for e->tau and mu->tau misidentification rates
  virtual double
  getSF_hadTauID_and_Iso(TauIDSFsys central_or_shift) const;

  virtual double
  getSF_eToTauFakeRate(FRet central_or_shift) const;

  virtual double
  getSF_muToTauFakeRate(FRmt central_or_shift) const;
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // data/MC corrections for jets to pass pileup jet ID
  double
  getSF_pileupJetID(pileupJetIDSFsys central_or_shift) const;
  //-----------------------------------------------------------------------------

protected:
  double
  getSF_leptonID_and_Iso(std::size_t numLeptons,
                         const std::vector<double> & lepton_pt,
                         const std::vector<double> & lepton_eta,
                         const std::vector<bool> & lepton_isGenMatched,
                         const std::vector<bool> & lepton_isTight,
                         bool sfForTightSelection,
                         const std::vector<lutWrapperBase *> & corrections,
                         int error_shift,
                         double recompSF) const;

  void
  initAntiEle_tauIDSFs(const std::string & era_str);

  void
  initAntiMu_tauIDSFs(const std::string & era_str);

  void
  init_tauIDSFs(const std::string & era_str,
                std::map<int, TauIDSFTool *> & tauIDSF_map,
                const std::string & tauID_str,
                int nof_levels);

  bool
  check_triggerSFsys_opt(TriggerSFsys central_or_shift) const;

  double
  comp_triggerSFsys_opt(double sf,
                        double sfErr,
                        TriggerSFsys central_or_shift) const;

  //-----------------------------------------------------------------------------
  // data/MC corrections for electron and muon identification and isolation efficiency,
  // including the cut on the ttH multilepton MVA

  // loose electron selection (RecoElectronSelectorLoose)
  std::vector<lutWrapperBase *> sfElectronID_and_Iso_loose_;
  // tight electron selection used in all channels except 2lss_1tau (RecoElectronSelectorTight with tightCharge_cut disabled)
  std::vector<lutWrapperBase *> sfElectronID_and_Iso_tight_to_loose_woTightCharge_;
  // tight electron selection specific to 2lss_1tau channel (RecoElectronSelectorTight with tightCharge_cut enabled)
  std::vector<lutWrapperBase *> sfElectronID_and_Iso_tight_to_loose_wTightCharge_;
  // errors for the above
  std::vector<lutWrapperBase *> sfElectronID_and_Iso_tight_to_loose_errors_up_;
  std::vector<lutWrapperBase *> sfElectronID_and_Iso_tight_to_loose_errors_down_;

  // loose muon selection (RecoMuonSelectorLoose)
  std::vector<lutWrapperBase *> sfMuonID_and_Iso_loose_;
  // tight muon selection used in all channels except 2lss_1tau (RecoMuonSelectorTight with tightCharge_cut disabled)
  std::vector<lutWrapperBase *> sfMuonID_and_Iso_tight_to_loose_woTightCharge_;
  // tight muon selection specific to 2lss_1tau channel (RecoMuonSelectorTight with tightCharge_cut enabled)
  std::vector<lutWrapperBase *> sfMuonID_and_Iso_tight_to_loose_wTightCharge_;
  // errors for the above
  std::vector<lutWrapperBase *> sfMuonID_and_Iso_tight_to_loose_errors_up_;
  std::vector<lutWrapperBase *> sfMuonID_and_Iso_tight_to_loose_errors_down_;
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // data/MC corrections for pileup jet ID
  lutWrapperBase * effPileupJetID_;
  lutWrapperBase * sfPileupJetID_eff_;
  lutWrapperBase * sfPileupJetID_eff_errors_;
  lutWrapperBase * mistagPileupJetID_;
  lutWrapperBase * sfPileupJetID_mistag_;
  lutWrapperBase * sfPileupJetID_mistag_errors_;
  //-----------------------------------------------------------------------------

  Era era_;

  std::map<std::string, TFile *> inputFiles_;

  int hadTauSelection_;
  TauID hadTauId_;
  std::string tauIDSF_str_;
  std::string tauIDSF_level_str_;

  int hadTauSelection_antiElectron_[4];
  int hadTauSelection_antiMuon_[4];

  TauIDSFTool * tauIdSFs_;
  std::map<int, TauIDSFTool *> tauIDSFs_antiEle_;
  std::map<int, TauIDSFTool *> tauIDSFs_antiMu_;
  bool applyHadTauSF_;
  bool isDEBUG_;

  pileupJetID pileupJetId_;

  bool recompTightSF_;
  double recompTightSF_el_woTightCharge_;
  double recompTightSF_mu_woTightCharge_;
  double recompTightSF_el_wTightCharge_;
  double recompTightSF_mu_wTightCharge_;

  std::size_t numLeptons_;
  std::vector<int> lepton_type_;
  std::vector<double> lepton_pt_;
  std::vector<double> lepton_cone_pt_;
  std::vector<double> lepton_eta_;
  std::size_t numElectrons_;
  std::vector<double> electron_pt_;
  std::vector<double> electron_cone_pt_;
  std::vector<double> electron_eta_;
  std::vector<bool> electron_isGenMatched_;
  std::vector<bool> electron_isTight_;
  std::size_t numMuons_;
  std::vector<double> muon_pt_;
  std::vector<double> muon_cone_pt_;
  std::vector<double> muon_eta_;
  std::vector<bool> muon_isGenMatched_;
  std::vector<bool> muon_isTight_;
  std::size_t numHadTaus_;
  std::vector<int> hadTau_genPdgId_;
  std::vector<double> hadTau_pt_;
  std::vector<double> hadTau_absEta_;
  std::size_t numJets_;
  std::vector<double> jet_pt_;
  std::vector<double> jet_eta_;
  std::vector<bool> jet_isPileup_;
  std::vector<bool> jet_passesPileupJetId_;
};

#endif // tthAnalysis_HiggsToTauTau_data_to_MC_corrections_Base_h
