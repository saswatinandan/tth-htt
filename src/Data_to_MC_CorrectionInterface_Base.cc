#include "tthAnalysis/HiggsToTauTau/interface/Data_to_MC_CorrectionInterface_Base.h"

#include "tthAnalysis/HiggsToTauTau/interface/lutAuxFunctions.h" // get_from_lut()
#include "tthAnalysis/HiggsToTauTau/interface/leptonTypes.h" // kElectron, kMuon
#include "tthAnalysis/HiggsToTauTau/interface/data_to_MC_corrections_auxFunctions.h" // aux::
#include "tthAnalysis/HiggsToTauTau/interface/sysUncertOptions.h" // pileupJetIDSFsys
#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h" // cmsException()
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // get_tau_id_wp_int()
#include "tthAnalysis/HiggsToTauTau/interface/RecoMuon.h" // RecoMuon
#include "tthAnalysis/HiggsToTauTau/interface/RecoElectron.h" // RecoElectron

#include "TauPOG/TauIDSFs/interface/TauIDSFTool.h" // TauIDSFTool

#include <TString.h> // Form()
#include <TFile.h> // TFile

#include <boost/algorithm/string/predicate.hpp> // boost::ends_with(), boost::starts_with()

#include <cassert> // assert()

namespace
{
  constexpr double
  square(double x)
  {
    return x * x;
  }
}

Data_to_MC_CorrectionInterface_Base::Data_to_MC_CorrectionInterface_Base(Era era, const edm::ParameterSet & cfg)
  : effPileupJetID_(nullptr)
  , sfPileupJetID_eff_(nullptr)
  , sfPileupJetID_eff_errors_(nullptr)
  , mistagPileupJetID_(nullptr)
  , sfPileupJetID_mistag_(nullptr)
  , sfPileupJetID_mistag_errors_(nullptr)
  , era_(era)
  , hadTauSelection_(-1)
  , hadTauId_(TauID::DeepTau2017v2VSjet)
  , tauIdSFs_(nullptr)
  , applyHadTauSF_(true)
  , isDEBUG_(cfg.exists("isDEBUG") ? cfg.getParameter<bool>("isDEBUG") : false)
  , pileupJetId_(kPileupJetID_disabled)
  , recompTightSF_(cfg.exists("lep_mva_wp") && cfg.getParameter<std::string>("lep_mva_wp") == "hh_multilepton")
  , recompTightSF_el_woTightCharge_(1.)
  , recompTightSF_mu_woTightCharge_(1.)
  , recompTightSF_el_wTightCharge_(1.)
  , recompTightSF_mu_wTightCharge_(1.)
  , numLeptons_(0)
  , numElectrons_(0)
  , numMuons_(0)
  , numHadTaus_(0)
  , numJets_(0)
{
  const std::string hadTauSelection_string = cfg.getParameter<std::string>("hadTauSelection");
  applyHadTauSF_ = hadTauSelection_string != "disabled";
  if(applyHadTauSF_)
  {
    setHadTauSelection(hadTauSelection_string);
  }

  for(int idxHadTau = 0; idxHadTau < 4; ++idxHadTau)
  {
    hadTauSelection_antiElectron_[idxHadTau] = -1;
    if(cfg.exists("hadTauSelection_antiElectron"))
    {
      hadTauSelection_antiElectron_[idxHadTau] = cfg.getParameter<int>("hadTauSelection_antiElectron");
    }
    else
    {
      const std::string cfgParName = "hadTauSelection_antiElectron" + aux::getHadTauIdxLabel(idxHadTau);
      if(cfg.exists(cfgParName))
      {
        hadTauSelection_antiElectron_[idxHadTau] = cfg.getParameter<int>(cfgParName);
      }
    }

    hadTauSelection_antiMuon_[idxHadTau] = -1;
    if(cfg.exists("hadTauSelection_antiMuon"))
    {
      hadTauSelection_antiMuon_[idxHadTau] = cfg.getParameter<int>("hadTauSelection_antiMuon");
    }
    else
    {
      const std::string cfgParName = "hadTauSelection_antiMuon" + aux::getHadTauIdxLabel(idxHadTau);
      if(cfg.exists(cfgParName))
      {
        hadTauSelection_antiMuon_[idxHadTau] = cfg.getParameter<int>(cfgParName);
      }
    }
  }

  if(applyHadTauSF_)
  {
    if(hadTauId_ == TauID::DeepTau2017v2VSjet)
    {
      tauIDSF_str_ = "DeepTau2017v2p1VSjet";
      tauIDSF_level_str_ = TauID_level_strings.at(TauID_levels.at(hadTauId_)).at(hadTauSelection_);
    }
    else if(hadTauId_ == TauID::MVAoldDM2017v2     ||
            hadTauId_ == TauID::MVAoldDMdR032017v2  )
    {
      tauIDSF_str_ = "MVAoldDM2017v2";
      tauIDSF_level_str_ = TauID_level_strings.at(TauID_levels.at(hadTauId_)).at(std::max(1, hadTauSelection_));
    }
    else
    {
      throw cmsException(this, __func__, __LINE__) << "Invalid tau ID: " << as_integer(hadTauId_);
    }
  }

  //-----------------------------------------------------------------------------
  // Efficiencies and mistag rates for pileup jet ID, provided by JetMET POG
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
  if ( cfg.exists("pileupJetID") )
  {
    pileupJetId_ = get_pileupJetID(cfg.getParameter<std::string>("pileupJetID"));
    if ( pileupJetId_ != kPileupJetID_disabled )
    {
      const std::string era_string = get_era(era_);
      std::string wp_string;
      // The encoding of the pileup jet ID working points is:
      //   puId==0 means 000: fail all PU ID;
      //   puId==4 means 100: pass loose ID, fail medium, fail tight;
      //   puId==6 means 110: pass loose and medium ID, fail tight;
      //   puId==7 means 111: pass loose, medium, tight ID. 
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID#miniAOD_and_nanoAOD
      if ( pileupJetId_ == kPileupJetID_loose )
      {
        wp_string = "L";
      }
      else if ( pileupJetId_ == kPileupJetID_medium )
      {
        wp_string = "M";
      }
      else if ( pileupJetId_ == kPileupJetID_tight )
      {
        wp_string = "T";
      }
      else
      {
        throw cmsException(this, __func__, __LINE__) << "Option 'pileupJetId_' set to an invalid value: " << pileupJetId_;
      }
      effPileupJetID_ = new lutWrapperTH2(
        inputFiles_,
        "tthAnalysis/HiggsToTauTau/data/pileupJetIdSF/effcyPUID_81Xtraining.root",
        Form("h2_eff_mc%s_%s", era_string.data(), wp_string.data()),
        lut::kXptYeta, 15., 50., lut::kLimit_and_Cut, -5.0, +5.0, lut::kLimit
      );
      sfPileupJetID_eff_ = new lutWrapperTH2(
        inputFiles_,
        "tthAnalysis/HiggsToTauTau/data/pileupJetIdSF/scalefactorsPUID_81Xtraining.root",
        Form("h2_eff_sf%s_%s", era_string.data(), wp_string.data()),
        lut::kXptYeta, 15., 50., lut::kLimit_and_Cut, -5.0, +5.0, lut::kLimit
      );
      sfPileupJetID_eff_errors_ = new lutWrapperTH2(
        inputFiles_,
        "tthAnalysis/HiggsToTauTau/data/pileupJetIdSF/scalefactorsPUID_81Xtraining.root",
        Form("h2_eff_sf%s_%s_Systuncty", era_string.data(), wp_string.data()),
        lut::kXptYeta, 15., 50., lut::kLimit_and_Cut, -5.0, +5.0, lut::kLimit
      );
      mistagPileupJetID_ = new lutWrapperTH2(
        inputFiles_,
        "tthAnalysis/HiggsToTauTau/data/pileupJetIdSF/effcyPUID_81Xtraining.root",
        Form("h2_mistag_mc%s_%s", era_string.data(), wp_string.data()),
        lut::kXptYeta, 15., 50., lut::kLimit_and_Cut, -5.0, +5.0, lut::kLimit
      );
      sfPileupJetID_mistag_ = new lutWrapperTH2(
        inputFiles_,
        "tthAnalysis/HiggsToTauTau/data/pileupJetIdSF/scalefactorsPUID_81Xtraining.root",
        Form("h2_mistag_sf%s_%s", era_string.data(), wp_string.data()),
        lut::kXptYeta, 15., 50., lut::kLimit_and_Cut, -5.0, +5.0, lut::kLimit
      );
      sfPileupJetID_mistag_errors_ = new lutWrapperTH2(
        inputFiles_,
        "tthAnalysis/HiggsToTauTau/data/pileupJetIdSF/scalefactorsPUID_81Xtraining.root",
        Form("h2_mistag_sf%s_%s_Systuncty", era_string.data(), wp_string.data()),
        lut::kXptYeta, 15., 50., lut::kLimit_and_Cut, -5.0, +5.0, lut::kLimit
      );
    }
  }
  //-----------------------------------------------------------------------------
  if(recompTightSF_)
  {
    if(era_ == Era::k2016)
    {
      // 2016 woTightCharge e   AvgSF_Def_el_woTightCharge: 0.796 +- 0.009;   AvgSF_New_el_woTightCharge: 0.875 +- 0.044
      // 2016 woTightCharge mu  AvgSF_Def_mu_woTightCharge: 0.922 +- 0.008;   AvgSF_New_mu_woTightCharge: 1.028 +- 0.040
      recompTightSF_el_woTightCharge_ = (1. - 0.875) / (1. - 0.796);
      recompTightSF_mu_woTightCharge_ = (1. - 1.028) / (1. - 0.922);
      
      // 2016 wTightCharge e   AvgSF_Def_el_wTightCharge: 0.790 +- 0.009;   AvgSF_New_el_wTightCharge: 0.867 +- 0.045
      // 2016 wTightCharge mu  AvgSF_Def_mu_wTightCharge: 0.922 +- 0.008;   AvgSF_New_mu_wTightCharge: 1.026 +- 0.040
      recompTightSF_el_wTightCharge_ = (1. - 0.867) / (1. - 0.790); 
      recompTightSF_mu_wTightCharge_ = (1. - 1.026) / (1. - 0.922); 
    }
    else if(era_ == Era::k2017)
    {
      // see: https://indico.cern.ch/event/961689/contributions/4047547/attachments/2114588/3557570/HHTo4W_3l_Updates_20201002_LooseLeptonSFCorrection_1.pdf
      // 2017 woTightCharge e   AvgSF_Def_el_woTightCharge: 0.755 +- 0.008;   AvgSF_New_el_woTightCharge: 0.886 +- 0.040
      // 2017 woTightCharge mu  AvgSF_Def_mu_woTightCharge: 0.881 +- 0.008;   AvgSF_New_mu_woTightCharge: 0.986 +- 0.034
      recompTightSF_el_woTightCharge_ = (1. - 0.886) / (1. - 0.755);
      recompTightSF_mu_woTightCharge_ = (1. - 0.986) / (1. - 0.881);
      
      // 2017 wTightCharge e   AvgSF_Def_el_wTightCharge: 0.750 +- 0.008;   AvgSF_New_el_wTightCharge: 0.871 +- 0.040 
      // 2017 wTightCharge mu  AvgSF_Def_mu_wTightCharge: 0.882 +- 0.008;   AvgSF_New_mu_wTightCharge: 0.981 +- 0.034
      recompTightSF_el_wTightCharge_ = (1. - 0.871) / (1. - 0.750);
      recompTightSF_mu_wTightCharge_ = (1. - 0.981) / (1. - 0.882); 
    }
    else if(era_ == Era::k2018)
    {
      // 2018 woTightCharge e   AvgSF_Def_el_woTightCharge: 0.834 +- 0.007;   AvgSF_New_el_woTightCharge: 0.929 +- 0.036 
      // 2018 woTightCharge mu  AvgSF_Def_mu_woTightCharge: 0.915 +- 0.006;   AvgSF_New_mu_woTightCharge: 0.956 +- 0.028 
      recompTightSF_el_woTightCharge_ = (1. - 0.929) / (1. - 0.834); 
      recompTightSF_mu_woTightCharge_ = (1. - 0.956) / (1. - 0.915);

      // 2018 wTightCharge e   AvgSF_Def_el_wTightCharge: 0.832 +- 0.007;   AvgSF_New_el_wTightCharge: 0.918 +- 0.036
      // 2018 wTightCharge mu  AvgSF_Def_mu_wTightCharge: 0.915 +- 0.006;   AvgSF_New_mu_wTightCharge: 0.950 +- 0.028
      recompTightSF_el_wTightCharge_ = (1. - 0.918) / (1. - 0.832);
      recompTightSF_mu_wTightCharge_ = (1. - 0.950) / (1. - 0.915); 
    }
    else
    {
      throw cmsException(this, __func__, __LINE__) << "Invalid era: " << as_integer(era_);
    }
  }
}

Data_to_MC_CorrectionInterface_Base::~Data_to_MC_CorrectionInterface_Base()
{
  aux::clearCollection(sfElectronID_and_Iso_loose_);
  aux::clearCollection(sfElectronID_and_Iso_tight_to_loose_woTightCharge_);
  aux::clearCollection(sfElectronID_and_Iso_tight_to_loose_wTightCharge_);
  aux::clearCollection(sfMuonID_and_Iso_loose_);
  aux::clearCollection(sfMuonID_and_Iso_tight_to_loose_woTightCharge_);
  aux::clearCollection(sfMuonID_and_Iso_tight_to_loose_wTightCharge_);
  for(auto & kv: inputFiles_)
  {
    delete kv.second;
  }
  if(applyHadTauSF_)
  {
    delete tauIdSFs_;
    for(auto & kv: tauIDSFs_antiEle_)
    {
      delete kv.second;
    }
    for(auto & kv: tauIDSFs_antiMu_)
    {
      delete kv.second;
    }
  }
  delete effPileupJetID_;
  delete sfPileupJetID_eff_;
  delete sfPileupJetID_eff_errors_;
  delete mistagPileupJetID_;
  delete sfPileupJetID_mistag_;
  delete sfPileupJetID_mistag_errors_;
}

void
Data_to_MC_CorrectionInterface_Base::setHadTauSelection(const std::string & hadTauSelection)
{
  hadTauSelection_ = get_tau_id_wp_int(hadTauSelection);
  hadTauId_ = get_tau_id_enum(hadTauSelection);
}

void
Data_to_MC_CorrectionInterface_Base::setLeptons(const std::vector<const RecoLepton *> & leptons,
                                                bool requireChargeMatch)
{
  numLeptons_ = 0;
  numMuons_ = 0;
  numElectrons_ = 0;

  lepton_type_.clear();
  lepton_pt_.clear();
  lepton_cone_pt_.clear();
  lepton_eta_.clear();

  muon_pt_.clear();
  muon_cone_pt_.clear();
  muon_eta_.clear();
  muon_isGenMatched_.clear();
  muon_isTight_.clear();
  
  electron_pt_.clear();
  electron_cone_pt_.clear();
  electron_eta_.clear();
  electron_isGenMatched_.clear();
  electron_isTight_.clear();
 
  for(const RecoLepton * const lepton: leptons)
  {
    lepton_pt_.push_back(lepton->pt());
    lepton_eta_.push_back(lepton->eta());
    lepton_cone_pt_.push_back(lepton->cone_pt());
    ++numLeptons_;

    const RecoMuon * const muon = dynamic_cast<const RecoMuon * const>(lepton);
    const RecoElectron * const electron = dynamic_cast<const RecoElectron * const>(lepton);
    if(muon)
    {
      lepton_type_.push_back(kMuon);

      muon_pt_.push_back(muon->pt());
      muon_cone_pt_.push_back(muon->cone_pt());
      muon_eta_.push_back(muon->eta());
      muon_isGenMatched_.push_back(muon->isGenMatched(requireChargeMatch));
      muon_isTight_.push_back(muon->isTight());
      ++numMuons_;
    }
    else if(electron)
    {
      lepton_type_.push_back(kElectron);

      electron_pt_.push_back(electron->pt());
      electron_cone_pt_.push_back(electron->cone_pt());
      electron_eta_.push_back(electron->eta());
      electron_isGenMatched_.push_back(electron->isGenMatched(requireChargeMatch) && electron->genLepton()); // [*]
      electron_isTight_.push_back(electron->isTight());
      ++numElectrons_;
      // [*] isGenMatched() can also return true if the reconstructed electron is matched to a generator level photon
    }
    else
    {
      assert(0);
    }
  }
}

void
Data_to_MC_CorrectionInterface_Base::setHadTaus(const std::vector<const RecoHadTau *> & hadTaus)
{
  numHadTaus_ = 0;
  hadTau_genPdgId_.clear();
  hadTau_pt_.clear();
  hadTau_absEta_.clear();
  for(const RecoHadTau * const hadTau: hadTaus)
  {
    hadTau_genPdgId_.push_back(getHadTau_genPdgId(hadTau));
    hadTau_pt_.push_back(hadTau->pt());
    hadTau_absEta_.push_back(hadTau->absEta());
    ++numHadTaus_;
  }
}

void
Data_to_MC_CorrectionInterface_Base::setJets(const std::vector<const RecoJet *> & jets)
{
  numJets_ = 0;
  jet_pt_.clear();
  jet_eta_.clear();
  jet_isPileup_.clear();
  jet_passesPileupJetId_.clear();
  for(const RecoJet * const jet: jets)
  {
    if(! jet->is_PUID_taggable())
    {
      continue;
    }
    jet_pt_.push_back(jet->pt());
    jet_eta_.push_back(jet->eta());
    jet_isPileup_.push_back(jet->genJet() ? false : true);
    jet_passesPileupJetId_.push_back(jet->passesPUID(pileupJetId_));
    ++numJets_;
  }
}


double
Data_to_MC_CorrectionInterface_Base::getSF_leptonTriggerEff(TriggerSFsys central_or_shift) const
{
  throw cmsException(this, __func__, __LINE__)
    << "Cannot call from base class"
  ;
}

double
Data_to_MC_CorrectionInterface_Base::getSF_leptonID_and_Iso(std::size_t numLeptons,
                                                            const std::vector<double> & lepton_pt,
                                                            const std::vector<double> & lepton_eta,
                                                            const std::vector<bool> & lepton_isGenMatched,
                                                            const std::vector<bool> & lepton_isTight,
                                                            bool sfForTightSelection,
                                                            const std::vector<lutWrapperBase *> & corrections,
                                                            int error_shift,
                                                            double recompSF) const
{
  if(isDEBUG_)
  {
    std::cout << get_human_line(this, __func__, __LINE__) << ": sfForTightSelection = " << sfForTightSelection << '\n';
  }
  double sf = 1.;
  for(std::size_t idxLepton = 0; idxLepton < numLeptons; ++idxLepton)
  {
    if(isDEBUG_)
    {
      std::cout << get_human_line(this, __func__, __LINE__)
                << "lepton_isGenMatched[idxLepton]: " << lepton_isGenMatched[idxLepton]
                << ", lepton_isTight[idxLepton]:" << lepton_isTight[idxLepton] << '\n'
      ;
    }
    if(! lepton_isGenMatched[idxLepton])
    {
      continue;
    }
    if(sfForTightSelection && (! lepton_isTight[idxLepton]))
    {
      continue;
    }
    const bool has_recompSF = std::fpclassify(recompSF) != FP_ZERO;
    const int error_shift_tmp = has_recompSF ? 0 : error_shift;
    const double sf_tmp = get_from_lut_err(corrections, lepton_pt[idxLepton], lepton_eta[idxLepton], error_shift_tmp, isDEBUG_);
    assert(sf_tmp > 0.);

    double corrFactor = 1.;
    if(has_recompSF)
    {
      // https://indico.cern.ch/event/961689/contributions/4047547/attachments/2114588/3557570/HHTo4W_3l_Updates_20201002_LooseLeptonSFCorrection_1.pdf
      assert(recompTightSF_);
      const double sf_recomp_tmp = 1. - (1. - sf_tmp) * recompSF;
      const double sf_recomp_shift = 0.5 * std::fabs(sf_tmp - sf_recomp_tmp);
      const double sf_recomp = sf_recomp_tmp + error_shift * sf_recomp_shift;
      corrFactor = sf_recomp / sf_tmp;
      if(isDEBUG_)
      {
        std::cout << get_human_line(this, __func__, __LINE__)
                  << "recompTightSF:: recompSF: " << recompSF
                  << ", sf: " << sf_tmp
                  << ", sf_corrected: " << sf_tmp * corrFactor << "\n";
      }
    }
    sf *= sf_tmp * corrFactor;
  }
  return sf;
}

void
Data_to_MC_CorrectionInterface_Base::initAntiEle_tauIDSFs(const std::string & era_str)
{
  return init_tauIDSFs(era_str, tauIDSFs_antiEle_, "antiEleMVA6", 5);
}

void
Data_to_MC_CorrectionInterface_Base::initAntiMu_tauIDSFs(const std::string & era_str)
{
  return init_tauIDSFs(era_str, tauIDSFs_antiMu_, "antiMu3", 2);
}

void
Data_to_MC_CorrectionInterface_Base::init_tauIDSFs(const std::string & era_str,
                                                   std::map<int, TauIDSFTool *> & tauIDSF_map,
                                                   const std::string & tauID_str,
                                                   int nof_levels)
{
  const std::vector<std::string> levels = TauID_level_strings.at(nof_levels);
  for(int level = 0; level < nof_levels; ++level)
  {
    tauIDSF_map[level + 1] = new TauIDSFTool(era_str, tauID_str, levels.at(level), false);
  }
}

bool
Data_to_MC_CorrectionInterface_Base::check_triggerSFsys_opt(TriggerSFsys central_or_shift) const
{
  if(central_or_shift == TriggerSFsys::central ||
     central_or_shift == TriggerSFsys::shiftUp ||
     central_or_shift == TriggerSFsys::shiftDown)
  {
    return true;
  }
  if(central_or_shift == TriggerSFsys::shift_2lssUp      ||
     central_or_shift == TriggerSFsys::shift_2lssDown    ||
     central_or_shift == TriggerSFsys::shift_2lssEEUp    ||
     central_or_shift == TriggerSFsys::shift_2lssEEDown  ||
     central_or_shift == TriggerSFsys::shift_2lssEMuUp   ||
     central_or_shift == TriggerSFsys::shift_2lssEMuDown ||
     central_or_shift == TriggerSFsys::shift_2lssMuMuUp  ||
     central_or_shift == TriggerSFsys::shift_2lssMuMuDown)
  {
    return numLeptons_ <= 2 && numHadTaus_ <= 2;
  }
  if(central_or_shift == TriggerSFsys::shift_3lUp ||
     central_or_shift == TriggerSFsys::shift_3lDown)
  {
    return numLeptons_ >= 3 && numHadTaus_ <= 1;
  }
  if(central_or_shift == TriggerSFsys::shift_1lEUp   ||
     central_or_shift == TriggerSFsys::shift_1lEDown ||
     central_or_shift == TriggerSFsys::shift_1lMuUp  ||
     central_or_shift == TriggerSFsys::shift_1lMuDown)
  {
    return numLeptons_ == 1 && numHadTaus_ == 0;
  }
  return false;
}

double
Data_to_MC_CorrectionInterface_Base::comp_triggerSFsys_opt(double sf,
                                                           double sfErr,
                                                           TriggerSFsys central_or_shift) const
{
  if(central_or_shift == TriggerSFsys::central ||
     ((central_or_shift == TriggerSFsys::shift_1lEUp      || central_or_shift == TriggerSFsys::shift_1lEDown     ) && (numElectrons_ != 1 || numMuons_ != 0)) ||
     ((central_or_shift == TriggerSFsys::shift_1lMuUp     || central_or_shift == TriggerSFsys::shift_1lMuDown    ) && (numElectrons_ != 0 || numMuons_ != 1)) ||
     ((central_or_shift == TriggerSFsys::shift_2lssEEUp   || central_or_shift == TriggerSFsys::shift_2lssEEDown  ) && (numElectrons_ != 2 || numMuons_ != 0)) ||
     ((central_or_shift == TriggerSFsys::shift_2lssEMuUp  || central_or_shift == TriggerSFsys::shift_2lssEMuDown ) && (numElectrons_ != 1 || numMuons_ != 1)) ||
     ((central_or_shift == TriggerSFsys::shift_2lssMuMuUp || central_or_shift == TriggerSFsys::shift_2lssMuMuDown) && (numElectrons_ != 0 || numMuons_ != 2)))
  {
    return sf;
  }
  else if(central_or_shift == TriggerSFsys::shiftUp          ||
          central_or_shift == TriggerSFsys::shift_1lEUp      ||
          central_or_shift == TriggerSFsys::shift_1lMuUp     ||
          central_or_shift == TriggerSFsys::shift_2lssUp     ||
          central_or_shift == TriggerSFsys::shift_2lssEEUp   ||
          central_or_shift == TriggerSFsys::shift_2lssEMuUp  ||
          central_or_shift == TriggerSFsys::shift_2lssMuMuUp ||
          central_or_shift == TriggerSFsys::shift_3lUp)
  {
    return sf * (1. + sfErr);
  }
  else if(central_or_shift == TriggerSFsys::shiftDown          ||
          central_or_shift == TriggerSFsys::shift_1lEDown      ||
          central_or_shift == TriggerSFsys::shift_1lMuDown     ||
          central_or_shift == TriggerSFsys::shift_2lssDown     ||
          central_or_shift == TriggerSFsys::shift_2lssEEDown   ||
          central_or_shift == TriggerSFsys::shift_2lssEMuDown  ||
          central_or_shift == TriggerSFsys::shift_2lssMuMuDown ||
          central_or_shift == TriggerSFsys::shift_3lDown)
  {
    return sf * (1. - sfErr);
  }
  throw cmsException(this, __func__, __LINE__) << "Invalid option: " << static_cast<int>(central_or_shift);
}

double
Data_to_MC_CorrectionInterface_Base::getSF_leptonID_and_Iso_loose(LeptonIDSFsys central_or_shift) const
{
  const bool sfForTightSelection = false;
  const double recompSF = 0.;
  if(isDEBUG_)
  {
    std::cout << get_human_line(this, __func__, __LINE__) << "Computing SF for electrons\n"
                "sfForTightSelection: " << sfForTightSelection << '\n'
    ;
  }
  int sf_el_error = 0;
  if(central_or_shift == LeptonIDSFsys::elLooseUp)
  {
    sf_el_error = +1;
  }
  else if(central_or_shift == LeptonIDSFsys::elLooseDown)
  {
    sf_el_error = -1;
  }
  const double sf_el = getSF_leptonID_and_Iso(
    numElectrons_, electron_pt_, electron_eta_, electron_isGenMatched_, electron_isTight_, sfForTightSelection, sfElectronID_and_Iso_loose_, sf_el_error, recompSF
  );
  if(isDEBUG_)
  {
    std::cout
      << "Final SF for electrons: " << sf_el << '\n'
      << get_human_line(this, __func__, __LINE__) << "Computing SF for muons\n"
    ;
  }

  int sf_mu_error = 0;
  if(central_or_shift == LeptonIDSFsys::muLooseUp)
  {
    sf_mu_error = +1;
  }
  else if(central_or_shift == LeptonIDSFsys::muLooseDown)
  {
    sf_mu_error = -1;
  }
  const double sf_mu = getSF_leptonID_and_Iso(
    numMuons_, muon_pt_, muon_eta_, muon_isGenMatched_, muon_isTight_, sfForTightSelection, sfMuonID_and_Iso_loose_, sf_mu_error, recompSF
  );
  const double sf = sf_el * sf_mu;
  if(isDEBUG_)
  {
    std::cout
      << "Final SF for muons: " << sf_mu << "\n"
         "Final SF: " << sf << '\n'
    ;
  }
  return sf;
}

double
Data_to_MC_CorrectionInterface_Base::getSF_leptonID_and_Iso_looseToFakeable() const
{
  return 1.; // CV: no data/MC corrections for "fakeable" leptons determined yet
}

double
Data_to_MC_CorrectionInterface_Base::getSF_leptonID_and_Iso_tight_to_loose_woTightCharge(LeptonIDSFsys central_or_shift) const
{
  const bool sfForTightSelection = true;
  const double recompSF = 0.;
  const int error_shift = 0;
  if(isDEBUG_)
  {
    std::cout << get_human_line(this, __func__, __LINE__) << "Computing SF for electrons\n"
                "sfForTightSelection: " << sfForTightSelection << '\n'
    ;
  }
  const double recompSF_el = recompTightSF_ ? recompTightSF_el_woTightCharge_ : 0.;
  int error_shift_el = 0;
  if(recompTightSF_)
  {
    if(central_or_shift == LeptonIDSFsys::elTightRecompUp)
    {
      error_shift_el = +1;
    }
    else if(central_or_shift == LeptonIDSFsys::elTightRecompDown)
    {
      error_shift_el = -1;
    }
  }
  double sf_el = getSF_leptonID_and_Iso(
    numElectrons_, electron_pt_, electron_eta_, electron_isGenMatched_, electron_isTight_, sfForTightSelection, sfElectronID_and_Iso_tight_to_loose_woTightCharge_, error_shift_el, recompSF_el
  );
  if(isDEBUG_)
  {
    std::cout << get_human_line(this, __func__, __LINE__) << "Computing SF errors for electrons\n";
  }
  if(central_or_shift == LeptonIDSFsys::elTightUp)
  {
    sf_el *= getSF_leptonID_and_Iso(
      numElectrons_, electron_pt_, electron_eta_, electron_isGenMatched_, electron_isTight_, sfForTightSelection, sfElectronID_and_Iso_tight_to_loose_errors_up_, error_shift, recompSF
    );
  }
  else if(central_or_shift == LeptonIDSFsys::elTightDown)
  {
    sf_el *= getSF_leptonID_and_Iso(
      numElectrons_, electron_pt_, electron_eta_, electron_isGenMatched_, electron_isTight_, sfForTightSelection, sfElectronID_and_Iso_tight_to_loose_errors_down_, error_shift, recompSF
    );
  }
  if(isDEBUG_)
  {
    std::cout
      << "Final SF for electrons: " << sf_el << '\n'
      << get_human_line(this, __func__, __LINE__) << "Computing SF for muons\n"
    ;
  }

  const double recompSF_mu = recompTightSF_ ? recompTightSF_mu_woTightCharge_ : 0.;
  int error_shift_mu = 0;
  if(recompTightSF_)
  {
    if(central_or_shift == LeptonIDSFsys::muTightRecompUp)
    {
      error_shift_mu = +1;
    }
    else if(central_or_shift == LeptonIDSFsys::muTightRecompDown)
    {
      error_shift_mu = -1;
    }
  }
  double sf_mu = getSF_leptonID_and_Iso(
    numMuons_, muon_pt_, muon_eta_, muon_isGenMatched_, muon_isTight_, sfForTightSelection, sfMuonID_and_Iso_tight_to_loose_woTightCharge_, error_shift_mu, recompSF_mu
  );
  if(isDEBUG_)
  {
    std::cout << get_human_line(this, __func__, __LINE__) << "Computing SF errors for muons\n";
  }
  if(central_or_shift == LeptonIDSFsys::muTightUp)
  {
    sf_mu *= getSF_leptonID_and_Iso(
      numMuons_, muon_pt_, muon_eta_, muon_isGenMatched_, muon_isTight_, sfForTightSelection, sfMuonID_and_Iso_tight_to_loose_errors_up_, error_shift, recompSF
    );
  }
  else if(central_or_shift == LeptonIDSFsys::muTightDown)
  {
    sf_mu *= getSF_leptonID_and_Iso(
      numMuons_, muon_pt_, muon_eta_, muon_isGenMatched_, muon_isTight_, sfForTightSelection, sfMuonID_and_Iso_tight_to_loose_errors_down_, error_shift, recompSF
    );
  }
  const double sf = sf_el * sf_mu;
  if(isDEBUG_)
  {
    std::cout
      << "Final SF for muons: " << sf_mu << "\n"
         "Final SF: " << sf << '\n'
    ;
  }
  return sf;
}

double
Data_to_MC_CorrectionInterface_Base::getSF_leptonID_and_Iso_tight_to_loose_wTightCharge(LeptonIDSFsys central_or_shift) const
{
  const bool sfForTightSelection = true;
  const double recompSF = 0.;
  const int error_shift = 0;
  if(isDEBUG_)
  {
    std::cout << get_human_line(this, __func__, __LINE__) << "Computing SF for electrons\n";
  }
  const double recompSF_el = recompTightSF_ ? recompTightSF_el_wTightCharge_ : 0.;
  int error_shift_el = 0;
  if(recompTightSF_)
  {
    if(central_or_shift == LeptonIDSFsys::elTightRecompUp)
    {
      error_shift_el = +1;
    }
    else if(central_or_shift == LeptonIDSFsys::elTightRecompDown)
    {
      error_shift_el = -1;
    }
  }  
  double sf_el = getSF_leptonID_and_Iso(
    numElectrons_, electron_pt_, electron_eta_, electron_isGenMatched_, electron_isTight_, sfForTightSelection, sfElectronID_and_Iso_tight_to_loose_wTightCharge_, error_shift_el, recompSF_el
  );
  if(isDEBUG_)
  {
    std::cout << get_human_line(this, __func__, __LINE__) << "Computing SF errors for electrons\n";
  }
  if(central_or_shift == LeptonIDSFsys::elTightUp)
  {
    sf_el *= getSF_leptonID_and_Iso(
      numElectrons_, electron_pt_, electron_eta_, electron_isGenMatched_, electron_isTight_, sfForTightSelection, sfElectronID_and_Iso_tight_to_loose_errors_up_, error_shift, recompSF
    );
  }
  else if(central_or_shift == LeptonIDSFsys::elTightDown)
  {
    sf_el *= getSF_leptonID_and_Iso(
      numElectrons_, electron_pt_, electron_eta_, electron_isGenMatched_, electron_isTight_, sfForTightSelection, sfElectronID_and_Iso_tight_to_loose_errors_down_, error_shift, recompSF
    );
  }
  if(isDEBUG_)
  {
    std::cout
      << "Final SF for electrons: " << sf_el << '\n'
      << get_human_line(this, __func__, __LINE__) << "Computing SF for muons\n"
    ;
  }

  const double recompSF_mu = recompTightSF_ ? recompTightSF_mu_wTightCharge_ : 0.;
  int error_shift_mu = 0;
  if(recompTightSF_)
  {
    if(central_or_shift == LeptonIDSFsys::muTightRecompUp)
    {
      error_shift_mu = +1;
    }
    else if(central_or_shift == LeptonIDSFsys::muTightRecompDown)
    {
      error_shift_mu = -1;
    }
  }
  double sf_mu = getSF_leptonID_and_Iso(
    numMuons_, muon_pt_, muon_eta_, muon_isGenMatched_, muon_isTight_, sfForTightSelection, sfMuonID_and_Iso_tight_to_loose_wTightCharge_, error_shift_mu, recompSF_mu
  );
  if(isDEBUG_)
  {
    std::cout << get_human_line(this, __func__, __LINE__) << "Computing SF errors for muons\n";
  }
  if(central_or_shift == LeptonIDSFsys::muTightUp)
  {
    sf_mu *= getSF_leptonID_and_Iso(
      numMuons_, muon_pt_, muon_eta_, muon_isGenMatched_, muon_isTight_, sfForTightSelection, sfMuonID_and_Iso_tight_to_loose_errors_up_, error_shift, recompSF
    );
  }
  else if(central_or_shift == LeptonIDSFsys::muTightDown)
  {
    sf_mu *= getSF_leptonID_and_Iso(
      numMuons_, muon_pt_, muon_eta_, muon_isGenMatched_, muon_isTight_, sfForTightSelection, sfMuonID_and_Iso_tight_to_loose_errors_down_, error_shift, recompSF
    );
  }
  const double sf = sf_el * sf_mu;
  if(isDEBUG_)
  {
    std::cout
      << "Final SF for muons: " << sf_mu << "\n"
         "Final SF: " << sf << '\n'
    ;
  }
  return sf;
}

double
Data_to_MC_CorrectionInterface_Base::getSF_hadTauID_and_Iso(TauIDSFsys central_or_shift) const
{
  double sf = 1.;
  if(applyHadTauSF_)
  {
    for(std::size_t idxHadTau = 0; idxHadTau < numHadTaus_; ++idxHadTau)
    {
      if(hadTau_genPdgId_[idxHadTau] == 15)
      {
        // because we use looser anti-e/mu DeepTau ID WPs compared to the WPs used in the tau ID SF measurement,
        // we have to add additional 3% or 15% uncertainty to the tau ID SF depending on the tau pT, as explained here:
        // https://indico.cern.ch/event/880308/contributions/3708888/attachments/1972379/3281439/NewsRun2SFsRecommendation.pdf
        const double sf_central = tauIdSFs_->getSFvsPT(hadTau_pt_[idxHadTau]);
        const double sf_up      = tauIdSFs_->getSFvsPT(hadTau_pt_[idxHadTau], "Up");
        const double sf_down    = tauIdSFs_->getSFvsPT(hadTau_pt_[idxHadTau], "Down");
        const double sf_unc = hadTau_pt_[idxHadTau] > 100. ? 0.15 : 0.03;
        const double sf_unc_up   = std::sqrt(square(sf_up   / sf_central - 1.) + square(sf_unc));
        const double sf_unc_down = std::sqrt(square(sf_down / sf_central - 1.) + square(sf_unc));
        switch(central_or_shift)
        {
          case TauIDSFsys::central:   sf *= sf_central;                      break;
          case TauIDSFsys::shiftUp:   sf *= sf_central * (1. + sf_unc_up);   break;
          case TauIDSFsys::shiftDown: sf *= sf_central * (1. - sf_unc_down); break;
        }
      }
    }
  }
  return sf;
}

double
Data_to_MC_CorrectionInterface_Base::getSF_eToTauFakeRate(FRet central_or_shift) const
{
  double sf = 1.;
  if(applyHadTauSF_)
  {
    for(std::size_t idxHadTau = 0; idxHadTau < numHadTaus_; ++idxHadTau)
    {
      const int hadTauSelection_antiElectron = hadTauSelection_antiElectron_[idxHadTau];
      const int hadTau_genPdgId = std::abs(hadTau_genPdgId_[idxHadTau]);
      if(hadTauSelection_antiElectron > 0 && hadTau_genPdgId  == 11)
      {
        if(! tauIDSFs_antiEle_.count(hadTauSelection_antiElectron))
        {
          throw cmsException(this, __func__, __LINE__) << "Anti-e SFs not initalized for WP " << hadTauSelection_antiElectron;
        }
        const TauIDSFTool * const tauIDSF_antiEle = tauIDSFs_antiEle_.at(hadTauSelection_antiElectron);

        switch(central_or_shift)
        {
          case FRet::shiftUp:   sf *= tauIDSF_antiEle->getSFvsEta(hadTau_absEta_[idxHadTau], 1, "Up");   break;
          case FRet::shiftDown: sf *= tauIDSF_antiEle->getSFvsEta(hadTau_absEta_[idxHadTau], 1, "Down"); break;
          case FRet::central:   sf *= tauIDSF_antiEle->getSFvsEta(hadTau_absEta_[idxHadTau], 1);         break;
        }
      }
    }
  }
  return sf;
}

double
Data_to_MC_CorrectionInterface_Base::getSF_muToTauFakeRate(FRmt central_or_shift) const 
{
  double sf = 1.;
  if(applyHadTauSF_)
  {
    for(std::size_t idxHadTau = 0; idxHadTau < numHadTaus_; ++idxHadTau)
    {
      const int hadTauSelection_antiMuon = hadTauSelection_antiMuon_[idxHadTau];
      const int hadTau_genPdgId = std::abs(hadTau_genPdgId_[idxHadTau]);
      if(hadTauSelection_antiMuon > 0 && hadTau_genPdgId  == 13)
      {
        if(! tauIDSFs_antiMu_.count(hadTauSelection_antiMuon))
        {
          throw cmsException(this, __func__, __LINE__) << "Anti-mu SFs not initalized for WP " << hadTauSelection_antiMuon;
        }
        const TauIDSFTool * const tauIDSF_antiMu = tauIDSFs_antiMu_.at(hadTauSelection_antiMuon);

        switch(central_or_shift)
        {
          case FRmt::shiftUp:   sf *= tauIDSF_antiMu->getSFvsEta(hadTau_absEta_[idxHadTau], 2, "Up");   break;
          case FRmt::shiftDown: sf *= tauIDSF_antiMu->getSFvsEta(hadTau_absEta_[idxHadTau], 2, "Down"); break;
          case FRmt::central:   sf *= tauIDSF_antiMu->getSFvsEta(hadTau_absEta_[idxHadTau], 2);         break;
        }
      }
    }
  }
  return sf;
}

double
Data_to_MC_CorrectionInterface_Base::getSF_pileupJetID(pileupJetIDSFsys central_or_shift) const
{
  // CV: Compute SF for efficiencies and mistag rates for jets to pass the pileup jet ID, following the recipe provided by the JetMET POG
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
  double prob_mc = 1.;
  double prob_data = 1.;
  if(isDEBUG_)
  {
    std::cout << get_human_line(this, __func__, __LINE__) << " -> sys unc = " << as_integer(central_or_shift) << '\n';
  }
  for(std::size_t idxJet = 0; idxJet < numJets_; ++idxJet)
  {
    if(isDEBUG_)
    {
      std::cout
        << get_human_line(this, __func__, __LINE__)
        << '#' << idxJet << " jet (pT = " << jet_pt_[idxJet] << ", eta = " << jet_eta_[idxJet] << ", "
           "is PU? = " << jet_isPileup_[idxJet] << ", passes PU? " << jet_passesPileupJetId_[idxJet] << ") -> "
      ;
    }
    double eff = 0.;
    double sf = 0.;
    if ( ! jet_isPileup_[idxJet] )
    {
      // apply efficiency to jets originating from hard-scatter interaction (ie real jets)
      eff = effPileupJetID_->getSF(jet_pt_[idxJet], jet_eta_[idxJet]);
      sf = sfPileupJetID_eff_->getSF(jet_pt_[idxJet], jet_eta_[idxJet]);
      const double sfErr = sfPileupJetID_eff_errors_->getSF(jet_pt_[idxJet], jet_eta_[idxJet]);
      if(isDEBUG_)
      {
        std::cout << "eff = " << eff << ", sf = " << sf << " +/- " << sfErr;
      }
      if      ( central_or_shift == pileupJetIDSFsys::effUp   ) sf += sfErr;
      else if ( central_or_shift == pileupJetIDSFsys::effDown ) sf -= sfErr;
    }
    else
    {
      // apply mistag rate to pileup jets
      eff = mistagPileupJetID_->getSF(jet_pt_[idxJet], jet_eta_[idxJet]);
      sf = sfPileupJetID_mistag_->getSF(jet_pt_[idxJet], jet_eta_[idxJet]);
      const double sfErr = sfPileupJetID_mistag_errors_->getSF(jet_pt_[idxJet], jet_eta_[idxJet]);
      if(isDEBUG_)
      {
        std::cout << "mistag = " << eff << ", sf = " << sf << " +/- " << sfErr;
      }
      if      ( central_or_shift == pileupJetIDSFsys::mistagUp   ) sf += sfErr;
      else if ( central_or_shift == pileupJetIDSFsys::mistagDown ) sf -= sfErr;
    }
    // CV: limit SF to the range [0..5], following the recommendation given by the JetMET POG on the twiki
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
    sf = std::clamp(sf, 0., 5.);
    if ( jet_passesPileupJetId_[idxJet] )
    {
      prob_mc *= eff;
      prob_data *= (sf*eff);
    }
    else
    {
      prob_mc *= (1.0 - eff);
      prob_data *= (1.0 - sf*eff);
    }
    if(isDEBUG_)
    {
      std::cout << " => sf = " << sf << ", prob_mc = " << prob_mc << ", prob_data = " << prob_data << '\n';
    }
  }
  const double sf = aux::compSF(prob_data, prob_mc);
  if(isDEBUG_)
  {
    std::cout
      << get_human_line(this, __func__, __LINE__)
      << "prob_mc = " << prob_mc << ", prob_data = " << prob_data << " => SF = " << sf << '\n'
    ;
  }
  return sf;
}
