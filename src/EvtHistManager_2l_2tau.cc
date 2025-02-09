#include "tthAnalysis/HiggsToTauTau/interface/EvtHistManager_2l_2tau.h"

#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // fillWithOverFlow(), fillWithOverFlow2d()
#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h" // cmsException()

EvtHistManager_2l_2tau::EvtHistManager_2l_2tau(const edm::ParameterSet & cfg)
  : HistManagerBase(cfg)
  , option_(kOption_undefined)
{
  const std::string option_string = cfg.getParameter<std::string>("option");
  if(option_string == "allHistograms")
  {
    option_ = kOption_allHistograms;
  }
  else if(option_string == "minimalHistograms")
  {
    option_ = kOption_minimalHistograms;
  }
  else
  {
    throw cmsException(__func__) << "Invalid Configuration parameter 'option' = " << option_string;
  }
  const std::vector<std::string> sysOpts_central = {
    "numElectrons",
    "numMuons",
    "numHadTaus",
    "numJets",
    "numBJets_loose",
    "numBJets_medium",
    "numBJets_loose_vs_numJets",
    "numBJets_medium_vs_numJets",
    "leptonPairCharge",
    "hadTauPairCharge",
    "mTauTauVis"
  };
  const std::vector<std::string> sysOpts_all = {
    "EventCounter",
    "mvaOutput_final"
  };
  for(const std::string & sysOpt: sysOpts_central)
  {
    central_or_shiftOptions_[sysOpt] = { "central" };
  }
  for(const std::string & sysOpt: sysOpts_all)
  {
    central_or_shiftOptions_[sysOpt] = { "*" };
  }
}

void
EvtHistManager_2l_2tau::bookHistograms(TFileDirectory & dir)
{
  if(option_ == kOption_allHistograms)
  {
    histogram_numElectrons_    = book1D(dir, "numElectrons",    "numElectrons",     5, -0.5,  +4.5);
    histogram_numMuons_        = book1D(dir, "numMuons",        "numMuons",         5, -0.5,  +4.5);
    histogram_numHadTaus_      = book1D(dir, "numHadTaus",      "numHadTaus",       5, -0.5,  +4.5);
    histogram_numJets_         = book1D(dir, "numJets",         "numJets",         20, -0.5, +19.5);
    histogram_numBJets_loose_  = book1D(dir, "numBJets_loose",  "numBJets_loose",  10, -0.5,  +9.5);
    histogram_numBJets_medium_ = book1D(dir, "numBJets_medium", "numBJets_medium", 10, -0.5,  +9.5);

    histogram_leptonPairCharge_ = book1D(dir, "leptonPairCharge", "leptonPairCharge", 5, -2.5, +2.5);
    histogram_hadTauPairCharge_ = book1D(dir, "hadTauPairCharge", "hadTauPairCharge", 5, -2.5, +2.5);
  }

  histogram_mTauTauVis_   = book1D(dir, "mTauTauVis",   "mTauTauVis",   40,  0.,  200.);
  histogram_EventCounter_ = book1D(dir, "EventCounter", "EventCounter",  1, -0.5,  +0.5);

  // X: binning by quantiles/era can be rounded to a unique one
  // 2018: [0.0, 0.7238368653359275, 1.0]
  // 2017: [0.0, 0.7085418998688564, 1.0]
  // 2016: [0.0, 0.6971249546108068, 1.0]
  Float_t binsx[3] = { 0.0, 0.255, 1.0 };
  histogram_final_ = book1D(dir, "mvaOutput_final",  "mvaOutput_final", 2, binsx);
}

void
EvtHistManager_2l_2tau::fillHistograms(const EvtHistManager_2l_2tau_Input & variables)
{
  const double evtWeightErr = 0.;
  const double & evtWeight = variables.evtWeight;

  if(option_ == kOption_allHistograms)
  {
    fillWithOverFlow(histogram_numElectrons_,    variables.numElectrons,    evtWeight, evtWeightErr);
    fillWithOverFlow(histogram_numMuons_,        variables.numMuons,        evtWeight, evtWeightErr);
    fillWithOverFlow(histogram_numHadTaus_,      variables.numHadTaus,      evtWeight, evtWeightErr);
    fillWithOverFlow(histogram_numJets_,         variables.numJets,         evtWeight, evtWeightErr);
    fillWithOverFlow(histogram_numBJets_loose_,  variables.numBJets_loose,  evtWeight, evtWeightErr);
    fillWithOverFlow(histogram_numBJets_medium_, variables.numBJets_medium, evtWeight, evtWeightErr);

    fillWithOverFlow(histogram_leptonPairCharge_, variables.leptonPairCharge, evtWeight, evtWeightErr);
    fillWithOverFlow(histogram_hadTauPairCharge_, variables.hadTauPairCharge, evtWeight, evtWeightErr);
  }

  fillWithOverFlow(histogram_mTauTauVis_,   variables.mTauTauVis, evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_EventCounter_, 0.,                   evtWeight, evtWeightErr);

  fillWithOverFlow(histogram_final_, variables.mvaOutput_legacy, evtWeight, evtWeightErr);
}
