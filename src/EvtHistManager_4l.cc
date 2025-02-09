#include "tthAnalysis/HiggsToTauTau/interface/EvtHistManager_4l.h"

#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // fillWithOverFlow(), getLogWeight()
#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h" // cmsException()

EvtHistManager_4l::EvtHistManager_4l(const edm::ParameterSet & cfg, bool isControlRegion)
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
    "numBJets_medium_vs_numJets"
  };
  std::vector<std::string> sysOpts_all = {
    "EventCounter"
  };
  if ( isControlRegion )
  {
    sysOpts_all.push_back("control");
  } else {
    sysOpts_all.push_back("mva_4l");
    sysOpts_all.push_back("mass_4L");
  }
  for(const std::string & sysOpt: sysOpts_central)
  {
    central_or_shiftOptions_[sysOpt] = { "central" };
  }
  for(const std::string & sysOpt: sysOpts_all)
  {
    central_or_shiftOptions_[sysOpt] = { "*" };
  }
}

const TH1 *
EvtHistManager_4l::getHistogram_EventCounter() const
{
  return histogram_EventCounter_;
}


void
EvtHistManager_4l::bookHistograms(TFileDirectory & dir)
{
  if(option_ == kOption_allHistograms)
  {
    histogram_numElectrons_    = book1D(dir, "numElectrons",    "numElectrons",     5, -0.5,  +4.5);
    histogram_numMuons_        = book1D(dir, "numMuons",        "numMuons",         5, -0.5,  +4.5);
    histogram_numJets_         = book1D(dir, "numJets",         "numJets",         20, -0.5, +19.5);
    histogram_numBJets_loose_  = book1D(dir, "numBJets_loose",  "numBJets_loose",  10, -0.5,  +9.5);
    histogram_numBJets_medium_ = book1D(dir, "numBJets_medium", "numBJets_medium", 10, -0.5,  +9.5);
  }

  histogram_EventCounter_ = book1D(dir, "EventCounter", "EventCounter", 1, -0.5, +0.5);
  Float_t bins_mass_4L[4] = { 70.,200.0,300.0,1000. };
  histogram_mass_4L_ = book1D(dir, "mass_4L", "mass_4L", 3, bins_mass_4L);
  Float_t binsx[3] = { 0.0, 0.55, 1.0 };
  histogram_mva_4l_ = book1D(dir, "mva_4l", "mva_4l", 2, binsx);

  histogram_ctrl_ = book1D(dir, "control", "control", 4, -0.5,  +3.5);
}

void
EvtHistManager_4l::fillHistograms(const EvtHistManager_4l_Input & variables)
{
  const double evtWeightErr = 0.;
  const double & evtWeight = variables.evtWeight;

  if(option_ == kOption_allHistograms)
  {
    fillWithOverFlow(histogram_numElectrons_,    variables.numElectrons,    evtWeight, evtWeightErr);
    fillWithOverFlow(histogram_numMuons_,        variables.numMuons,        evtWeight, evtWeightErr);
    fillWithOverFlow(histogram_numJets_,         variables.numJets,         evtWeight, evtWeightErr);
    fillWithOverFlow(histogram_numBJets_loose_,  variables.numBJets_loose,  evtWeight, evtWeightErr);
    fillWithOverFlow(histogram_numBJets_medium_, variables.numBJets_medium, evtWeight, evtWeightErr);
  }

  fillWithOverFlow(histogram_mass_4L_, variables.mass_4L,       evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_mva_4l_,  variables.mva_4l,        evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_ctrl_,    variables.ctrl_category, evtWeight, evtWeightErr);

  fillWithOverFlow(histogram_EventCounter_, 0., evtWeight, evtWeightErr);
}
