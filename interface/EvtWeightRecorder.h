#ifndef tthAnalysis_HiggsToTauTau_EvtWeightRecorder_h
#define tthAnalysis_HiggsToTauTau_EvtWeightRecorder_h

#include "tthAnalysis/HiggsToTauTau/interface/GenParticle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h" // edm::VParameterSet

// forward declarations
class L1PreFiringWeightReader;
class LHEInfoReader;
class PSWeightReader;
class EventInfo;
class RecoJet;
class Data_to_MC_CorrectionInterface_Base;
class JetToTauFakeRateInterface;
class LeptonFakeRateInterface;
class RecoLepton;
class RecoHadTau;
class EvtWeightManager;
class DYMCReweighting;
class DYMCNormScaleFactors;
class Data_to_MC_CorrectionInterface_0l_2tau_trigger;
class Data_to_MC_CorrectionInterface_1l_1tau_trigger;
class Data_to_MC_CorrectionInterface_1l_2tau_trigger;
class BtagSFRatioFacility;
class HHWeightInterfaceLO;
class HHWeightInterfaceNLO;
class LHEVpt_LOtoNLO;
class SubjetBtagSF;

enum class L1PreFiringWeightSys;
enum class PDFSys;
enum class PUsys;
enum class TriggerSFsys;
enum class LeptonIDSFsys;
enum class TauIDSFsys;
enum class FRet;
enum class FRmt;
enum class EWKJetSys;
enum class EWKBJetSys;
enum class pileupJetIDSFsys;
enum class LHEVptSys;
enum class SubjetBtagSys;

class EvtWeightRecorder
{
public:
  EvtWeightRecorder();
  EvtWeightRecorder(const std::vector<std::string> & central_or_shifts,
                    const std::string & central_or_shift,
                    bool isMC);
  virtual ~EvtWeightRecorder() {}

  double
  get(const std::string & central_or_shift,
      const std::string & bin = "") const;

  virtual double
  get_inclusive(const std::string & central_or_shift,
                const std::string & bin = "") const;

  double
  get_genWeight() const;

  double
  get_hhWeight() const;

  double
  get_hhWeight_lo() const;

  double
  get_hhWeight_nlo() const;

  double
  get_rescaling() const;

  double
  get_auxWeight(const std::string & central_or_shift) const;

  double
  get_lumiScale(const std::string & central_or_shift,
                const std::string & bin = "") const;

  double
  get_prescaleWeight() const; 
  
  double
  get_btagSFRatio(const std::string & central_or_shift) const;

  double
  get_nom_tH_weight(const std::string & central_or_shift) const;

  double
  get_puWeight(const std::string & central_or_shift) const;

  double
  get_pileupJetIDSF(const std::string & central_or_shift) const;

  double
  get_l1PreFiringWeight(const std::string & central_or_shift) const;

  double
  get_lheScaleWeight(const std::string & central_or_shift) const;

  double
  get_pdfWeight(const std::string & central_or_shift) const;

  double
  get_pdfMemberWeight(const std::string & central_or_shift) const;

  double
  get_psWeight(const std::string & central_or_shift) const;

  double
  get_leptonSF() const;

  double
  get_leptonIDSF_recoToLoose(const std::string & central_or_shift) const;

  double
  get_leptonIDSF_looseToTight(const std::string & central_or_shift) const;

  double
  get_leptonIDSF(const std::string & central_or_shift) const;

  double
  get_chargeMisIdProb() const;

  double
  get_dyBgrWeight() const; // CV: weight for data-driven Drell-Yan background estimation used in HH->bbWW dilepton channel

  double
  get_data_to_MC_correction(const std::string & central_or_shift) const;

  double
  get_btag(const std::string & central_or_shift) const;

  double
  get_ewk_jet(const std::string & central_or_shift) const;

  double
  get_ewk_bjet(const std::string & central_or_shift) const;

  double
  get_dy_rwgt(const std::string & central_or_shift) const;

  double
  get_dy_norm(const std::string & central_or_shift) const;

  double
  get_toppt_rwgt(const std::string & central_or_shift) const;

  double
  get_LHEVpt(const std::string & central_or_shift) const;

  double
  get_subjetBtagSF(const std::string & central_or_shift) const;
  
  double
  get_sf_triggerEff(const std::string & central_or_shift) const;

  virtual double
  get_tauSF(const std::string & central_or_shift) const;

  double
  get_FR(const std::string & central_or_shift) const;

  void
  record_genWeight(const EventInfo & eventInfo,
                   bool use_sign_only = false);

  void
  record_auxWeight(const EvtWeightManager * const evtWeightManager);

  void
  record_dy_rwgt(const DYMCReweighting * const dyReweighting,
                 const std::vector<GenParticle> & genTauLeptons);

  void
  record_dy_norm(const DYMCNormScaleFactors * const dyNormScaleFactors,
                 const std::vector<GenParticle> & genTauLeptons,
                 int nJets,
                 int nBLoose,
                 int nBMedium);

  void
  record_toppt_rwgt(double sf);
  
  void
  record_lumiScale(const edm::VParameterSet & lumiScales);

  void
  record_btagSFRatio(const BtagSFRatioFacility * const btagSFRatioFacility,
                     int nselJets);

  void
  record_rescaling(double rescaling);

  void
  record_nom_tH_weight(const EventInfo * const eventInfo);

  void
  record_leptonSF(double weight);

  void
  record_leptonIDSF_recoToLoose(const Data_to_MC_CorrectionInterface_Base * const dataToMCcorrectionInterface);

  void
  record_leptonIDSF_looseToTight(const Data_to_MC_CorrectionInterface_Base * const dataToMCcorrectionInterface,
                                 bool woTightCharge = true);

  void
  record_chargeMisIdProb(double weight);

  // CV: weight for data-driven Drell-Yan background estimation used in HH->bbWW dilepton channel
  void
  record_dyBgrWeight(double weight); 

  void
  record_prescale(double weight);

  // CV: weight specific to LO HH MC samples, computed by HHWeightInterfaceLO->getWeight("SM",...) method
  void
  record_hhWeight_lo(double weight);

  void
  record_hhWeight_lo(const HHWeightInterfaceLO * const HHWeightLO_calc,
                     const EventInfo & eventInfo,
                     bool isDEBUG = false);

  // CV: weight specific to LO HH MC samples, computed by HHWeightInterfaceNLO->getWeight("SM",...) method
  void
  record_hhWeight_nlo(double weight);

  void
  record_hhWeight_nlo(const HHWeightInterfaceNLO * const HHWeightNLO_calc,
                      const EventInfo & eventInfo,
                      bool isDEBUG = false);

  void
  record_l1PrefireWeight(const L1PreFiringWeightReader * const l1PreFiringWeightReader);

  void
  record_lheScaleWeight(const LHEInfoReader * const lheInfoReader);

  void
  record_pdfWeight(const LHEInfoReader * const lheInfoReader);

  void
  record_pdfMembers(const LHEInfoReader * const lheInfoReader,
                    const std::map<std::string, int> & pdf_map);

  void
  record_psWeight(const PSWeightReader * const psWeightReader);

  void
  record_puWeight(const EventInfo * const eventInfo);

  void
  record_pileupJetIDSF(const Data_to_MC_CorrectionInterface_Base * const dataToMCcorrectionInterface);

  void
  record_btagWeight(const std::vector<const RecoJet *> & jets);

  void
  record_ewk_jet(const std::vector<const RecoJet *> & jets);

  void
  record_ewk_bjet(const std::vector<const RecoJet *> & bjets);

  void
  record_LHEVpt(const LHEVpt_LOtoNLO * const lhe_vpt);

  void
  record_subjetBtagSF(SubjetBtagSF * const subjetBtagSF);

  void
  record_leptonTriggerEff(const Data_to_MC_CorrectionInterface_Base * const dataToMCcorrectionInterface);

  void
  record_tauTriggerEff(const Data_to_MC_CorrectionInterface_0l_2tau_trigger * const dataToMCcorrectionInterface_0l_2tau_trigger);

  void
  record_tauTriggerEff(const Data_to_MC_CorrectionInterface_1l_1tau_trigger * const dataToMCcorrectionInterface_1l_1tau_trigger);

  void
  record_tauTriggerEff(const Data_to_MC_CorrectionInterface_1l_2tau_trigger * const dataToMCcorrectionInterface_1l_2tau_trigger);

  void
  record_hadTauID_and_Iso(const Data_to_MC_CorrectionInterface_Base * const dataToMCcorrectionInterface);

  void
  record_eToTauFakeRate(const Data_to_MC_CorrectionInterface_Base * const dataToMCcorrectionInterface);

  void
  record_muToTauFakeRate(const Data_to_MC_CorrectionInterface_Base * const dataToMCcorrectionInterface);

  void
  record_jetToTau_FR_lead(const JetToTauFakeRateInterface * const jetToTauFakeRateInterface,
                          const RecoHadTau * const hadTau_lead);

  void
  record_jetToTau_FR_sublead(const JetToTauFakeRateInterface * const jetToTauFakeRateInterface,
                             const RecoHadTau * const hadTau_sublead);

  void
  record_jetToTau_SF_lead(const JetToTauFakeRateInterface * const jetToTauFakeRateInterface,
                          const RecoHadTau * const hadTau_lead);

  void
  record_jetToTau_SF_sublead(const JetToTauFakeRateInterface * const jetToTauFakeRateInterface,
                             const RecoHadTau * const hadTau_sublead);

  void
  record_jetToLepton_FR_lead(const LeptonFakeRateInterface * const leptonFakeRateInterface,
                             const RecoLepton * const lepton_lead);

  void
  record_jetToLepton_FR_sublead(const LeptonFakeRateInterface * const leptonFakeRateInterface,
                                const RecoLepton * const lepton_sublead);

  void
  record_jetToLepton_FR_third(const LeptonFakeRateInterface * const leptonFakeRateInterface,
                              const RecoLepton * const lepton_third);

  void
  record_jetToLepton_FR_fourth(const LeptonFakeRateInterface * const leptonFakeRateInterface,
                               const RecoLepton * const lepton_fourth);

  void
  compute_FR_2l2tau(bool passesTight_lepton_lead,
                    bool passesTight_lepton_sublead,
                    bool passesTight_hadTau_lead,
                    bool passesTight_hadTau_sublead);

  void
  compute_FR_2l1tau(bool passesTight_lepton_lead,
                    bool passesTight_lepton_sublead,
                    bool passesTight_hadTau);

  void
  compute_FR_2l(bool passesTight_lepton_lead,
                bool passesTight_lepton_sublead);

  void
  compute_FR_3l(bool passesTight_lepton_lead,
                bool passesTight_lepton_sublead,
                bool passesTight_lepton_third);

  void
  compute_FR_4l(bool passesTight_lepton_lead,
                bool passesTight_lepton_sublead,
                bool passesTight_lepton_third,
                bool passesTight_lepton_fourth);

  void
  compute_FR_3l1tau(bool passesTight_lepton_lead,
                    bool passesTight_lepton_sublead,
                    bool passesTight_lepton_third,
                    bool passesTight_hadTau);

  void
  compute_FR_1l2tau(bool passesTight_lepton,
                    bool passesTight_hadTau_lead,
                    bool passesTight_hadTau_sublead);

  void
  compute_FR_1l1tau(bool passesTight_lepton,
                    bool passesTight_hadTau);

  void
  compute_FR_1l(bool passesTight_lepton);

  void
  compute_FR_1tau(bool passesTight_hadTau);

  void
  compute_FR_2tau(bool passesTight_hadTau_lead,
                  bool passesTight_hadTau_sublead);

  double
  get_jetToLepton_FR_lead(const std::string & central_or_shift);

  double
  get_jetToLepton_FR_sublead(const std::string & central_or_shift);

  double
  get_jetToLepton_FR_third(const std::string & central_or_shift);

  double
  get_jetToLepton_FR_fourth(const std::string & central_or_shift);

  double
  get_jetToTau_FR_lead(const std::string & central_or_shift);

  double
  get_jetToTau_FR_sublead(const std::string & central_or_shift);

  friend std::ostream &
  operator<<(std::ostream & os,
             const EvtWeightRecorder & evtWeightRecorder);

protected:
  void
  record_jetToLepton_FR(const LeptonFakeRateInterface * const leptonFakeRateInterface,
                        const RecoLepton * const lepton,
                        std::map<int, double> & weights_FR_lepton);

  bool isMC_;
  double genWeight_;
  std::map<std::string, double> auxWeight_;
  std::map<std::string, std::map<std::string, double>> lumiScale_;
  std::map<std::string, double> nom_tH_weight_;
  std::map<std::string, double> btagSFRatio_;
  double leptonSF_;
  double chargeMisIdProb_;
  double dyBgrWeight_;
  double prescale_;
  double hhWeight_lo_;
  double hhWeight_nlo_;
  double rescaling_;
  std::string central_or_shift_;
  std::vector<std::string> central_or_shifts_;

  std::map<L1PreFiringWeightSys, double> weights_l1PreFiring_;
  std::map<int, double> weights_lheScale_;
  std::map<PDFSys, double> weights_pdf_;
  std::map<std::string, double> weights_pdf_members_;
  std::map<int, double> weights_partonShower_;
  std::map<PUsys, double> weights_pu_;
  std::map<int, double> weights_dy_norm_;
  std::map<int, double> weights_dy_rwgt_;
  std::map<int, double> weights_toppt_rwgt_;
  std::map<LHEVptSys, double> weights_lhe_vpt_;
  std::map<SubjetBtagSys, double> weights_subjet_btag_;
  std::map<int, double> weights_btag_;
  std::map<pileupJetIDSFsys, double> weights_puJetIDSF_;
  std::map<TriggerSFsys, double> weights_leptonTriggerEff_;
  std::map<TriggerSFsys, double> weights_tauTriggerEff_;
  std::map<LeptonIDSFsys, double> weights_leptonID_and_Iso_recoToLoose_;
  std::map<LeptonIDSFsys, double> weights_leptonID_and_Iso_looseToTight_;
  std::map<TauIDSFsys, double> weights_hadTauID_and_Iso_;
  std::map<FRet, double> weights_eToTauFakeRate_;
  std::map<FRmt, double> weights_muToTauFakeRate_;
  std::map<int, double> weights_FR_hadTau_lead_;
  std::map<int, double> weights_FR_hadTau_sublead_;
  std::map<int, double> weights_SF_hadTau_lead_;
  std::map<int, double> weights_SF_hadTau_sublead_;
  std::map<int, double> weights_FR_lepton_lead_;
  std::map<int, double> weights_FR_lepton_sublead_;
  std::map<int, double> weights_FR_lepton_third_;
  std::map<int, double> weights_FR_lepton_fourth_;
  std::map<std::string, double> weights_FR_;
  std::map<EWKJetSys, double> weights_ewk_jet_;
  std::map<EWKBJetSys, double> weights_ewk_bjet_;
};

#endif
