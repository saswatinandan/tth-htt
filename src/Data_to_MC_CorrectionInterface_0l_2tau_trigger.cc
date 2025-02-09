#include "tthAnalysis/HiggsToTauTau/interface/Data_to_MC_CorrectionInterface_0l_2tau_trigger.h"

#include "tthAnalysis/HiggsToTauTau/interface/leptonTypes.h" // kElectron, kMuon
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // Era::k*, get_era()
#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h" // cmsException()
#include "tthAnalysis/HiggsToTauTau/interface/data_to_MC_corrections_auxFunctions.h" // aux::
#include "tthAnalysis/HiggsToTauTau/interface/sysUncertOptions.h" // TriggerSFsys
#include "tthAnalysis/HiggsToTauTau/interface/lutAuxFunctions.h" // get_from_lut()

#include "TauAnalysisTools/TauTriggerSFs/interface/TauTriggerSFs2017.h" // TauTriggerSFs2017

#include <TString.h> // Form()
#include <TMath.h> // TMath::Abs()

#include <boost/algorithm/string/predicate.hpp> // boost::ends_with()

#include <cassert> // assert()

Data_to_MC_CorrectionInterface_0l_2tau_trigger::Data_to_MC_CorrectionInterface_0l_2tau_trigger(const edm::ParameterSet & cfg)
  : era_str_(cfg.getParameter<std::string>("era"))
  , era_(get_era(era_str_))
  , hadTauSelection_(cfg.getParameter<std::string>("hadTauSelection"))
  , isDEBUG_(cfg.exists("isDEBUG") ? cfg.getParameter<bool>("isDEBUG") : false)
  , allowedDecayModes_({ 0, 1, 2, 10, 11 })
  , hadTau1_genPdgId_(0)
  , hadTau1_pt_(0.)
  , hadTau1_eta_(0.)
  , hadTau1_phi_(0.)
  , hadTau1_decayMode_(0)
  , hadTau2_genPdgId_(0)
  , hadTau2_pt_(0.)
  , hadTau2_eta_(0.)
  , hadTau2_phi_(0.)
  , hadTau2_decayMode_(0)
  , effTrigger_tauLeg_(era_str_, hadTauSelection_, TauTriggerType::DiTau)
{}

Data_to_MC_CorrectionInterface_0l_2tau_trigger::~Data_to_MC_CorrectionInterface_0l_2tau_trigger()
{}

void
Data_to_MC_CorrectionInterface_0l_2tau_trigger::setTriggerBits(bool isTriggered_2tau)
{
  isTriggered_2tau_ = isTriggered_2tau;
}

void
Data_to_MC_CorrectionInterface_0l_2tau_trigger::setHadTaus(const RecoHadTau * const hadTau1,
                                                           const RecoHadTau * const hadTau2)
{
  hadTau1_pt_ = hadTau1->pt();
  hadTau1_eta_ = hadTau1->eta();
  hadTau1_phi_ = hadTau1->phi();
  hadTau1_decayMode_ = hadTau1->decayMode();

  hadTau2_pt_ = hadTau2->pt();
  hadTau2_eta_ = hadTau2->eta();
  hadTau2_phi_ = hadTau2->phi();
  hadTau2_decayMode_ = hadTau2->decayMode();
}

double
Data_to_MC_CorrectionInterface_0l_2tau_trigger::getSF_triggerEff(TriggerSFsys central_or_shift) const
{
  if(isDEBUG_)
  {
    std::cout << get_human_line(this, __func__, __LINE__) << '\n';
  }
  assert(check_triggerSFsys_opt(central_or_shift));

  double eff_2tau_tauLeg1_data = 0.;
  double eff_2tau_tauLeg1_mc   = 0.;
  double eff_2tau_tauLeg2_data = 0.;
  double eff_2tau_tauLeg2_mc   = 0.;

  if(std::fabs(hadTau1_eta_) <= 2.1 && aux::hasDecayMode(allowedDecayModes_, hadTau1_decayMode_))
  {
    eff_2tau_tauLeg1_data = effTrigger_tauLeg_.getTauTriggerEvalData(central_or_shift, hadTau1_pt_, hadTau1_eta_, hadTau1_phi_, hadTau1_decayMode_);
    eff_2tau_tauLeg1_mc   = effTrigger_tauLeg_.getTauTriggerEvalMC  (central_or_shift, hadTau1_pt_, hadTau1_eta_, hadTau1_phi_, hadTau1_decayMode_);
  }

  if(std::fabs(hadTau2_eta_) <= 2.1 && aux::hasDecayMode(allowedDecayModes_, hadTau2_decayMode_))
  {
    eff_2tau_tauLeg2_data = effTrigger_tauLeg_.getTauTriggerEvalData(central_or_shift, hadTau2_pt_, hadTau2_eta_, hadTau2_phi_, hadTau2_decayMode_);
    eff_2tau_tauLeg2_mc   = effTrigger_tauLeg_.getTauTriggerEvalMC  (central_or_shift, hadTau2_pt_, hadTau2_eta_, hadTau2_phi_, hadTau2_decayMode_);
  }

  if(isDEBUG_)
  {
    std::cout << "hadTau (lead):    pT = " << hadTau1_pt_ << ", eta = " << hadTau1_eta_ << "\n"
                 "hadTau (sublead): pT = " << hadTau2_pt_ << ", eta = " << hadTau2_eta_ << "\n"
                 "eff (tau1): data = " << eff_2tau_tauLeg1_data << ", MC = " << eff_2tau_tauLeg1_mc << "\n"
                 "eff (tau2): data = " << eff_2tau_tauLeg2_data << ", MC = " << eff_2tau_tauLeg2_mc << '\n'
    ;
  }

  double sf = 1.;
  if(isTriggered_2tau_)
  {
    // case 2: ditau trigger fires

    const double eff_data = eff_2tau_tauLeg1_data * eff_2tau_tauLeg2_data;
    const double eff_mc   = eff_2tau_tauLeg1_mc   * eff_2tau_tauLeg2_mc;
    if ( eff_data < 1.e-3 && eff_mc < 1.e-3 )
    {
      std::cout << "tauLeg1: pT = " << hadTau1_pt_ << ", eta = " << hadTau1_eta_ << std::endl;
      std::cout << "tauLeg2: pT = " << hadTau2_pt_ << ", eta = " << hadTau2_eta_ << std::endl;
      std::cout << " eff_data = " << eff_data << ": eff_2tau_tauLeg1_data = " << eff_2tau_tauLeg1_data << ", eff_2tau_tauLeg1_data = " << eff_2tau_tauLeg2_data << std::endl;
      std::cout << " eff_mc = " << eff_mc << ": eff_2tau_tauLeg1_mc = " << eff_2tau_tauLeg1_mc << ", eff_2tau_tauLeg1_mc = " << eff_2tau_tauLeg2_mc << std::endl;
    }
    sf = aux::compSF(eff_data, eff_mc);
    if(isDEBUG_)
    {
      std::cout << "case 2: ditau trigger fires\n"
                   " eff: data = " << eff_data << ", MC = " << eff_mc << " --> SF = " << sf << '\n'
      ;
    }
  }
  else
  {
    // case 1: ditau trigger doesn't fire
    // (SF doesn't matter, as event does not pass event selection)

    sf = 0.;
    if(isDEBUG_)
    {
      std::cout << "case 1: ditau trigger doesn't fire\n"
                   "--> setting SF = " << sf << '\n'
      ;
    }
  }

  return sf;
}

bool
Data_to_MC_CorrectionInterface_0l_2tau_trigger::check_triggerSFsys_opt(TriggerSFsys central_or_shift) const
{
  return
    central_or_shift == TriggerSFsys::central          ||
    central_or_shift == TriggerSFsys::shiftUp          ||
    central_or_shift == TriggerSFsys::shiftDown        ||
    central_or_shift == TriggerSFsys::shift_0l2tauUp   ||
    central_or_shift == TriggerSFsys::shift_0l2tauDown
  ;
}
