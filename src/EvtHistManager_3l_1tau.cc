#include "tthAnalysis/HiggsToTauTau/interface/EvtHistManager_3l_1tau.h"

#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // fillWithOverFlow(), getLogWeight()
#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h" // cmsException()

EvtHistManager_3l_1tau::EvtHistManager_3l_1tau(const edm::ParameterSet & cfg)
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
    "mTauTauVis",
    "memOutput_isValid",
    "memOutput_errorFlag",
    "memOutput_logWeight_ttH",
    "memOutput_logWeight_ttZ",
    "memOutput_logWeight_ttH_hww",
    "memOutput_LR",
    "mem_logCPUTime",
    "mem_logRealTime",
  };
  const std::vector<std::string> sysOpts_all = {
    "mvaOutput_legacy",
    "EventCounter",
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

const TH1 *
EvtHistManager_3l_1tau::getHistogram_EventCounter() const
{
  return histogram_EventCounter_;
}

void
EvtHistManager_3l_1tau::bookHistograms(TFileDirectory & dir)
{
  if(option_ == kOption_allHistograms)
  {
    histogram_numElectrons_    = book1D(dir, "numElectrons",    "numElectrons",     5, -0.5,  +4.5);
    histogram_numMuons_        = book1D(dir, "numMuons",        "numMuons",         5, -0.5,  +4.5);
    histogram_numHadTaus_      = book1D(dir, "numHadTaus",      "numHadTaus",       5, -0.5,  +4.5);
    histogram_numJets_         = book1D(dir, "numJets",         "numJets",         20, -0.5, +19.5);
    histogram_numBJets_loose_  = book1D(dir, "numBJets_loose",  "numBJets_loose",  10, -0.5,  +9.5);
    histogram_numBJets_medium_ = book1D(dir, "numBJets_medium", "numBJets_medium", 10, -0.5,  +9.5);
  }

  // X: binning by quantiles/era -- we can round to the same for all eras
  //2016: [0.0, 0.36713849092205514, 0.4889022577306826, 0.5865429137829421, 1.0]
  //2017: [0.0, 0.36597589479860454, 0.47015027737107856, 0.5638581399225489, 1.0]
  //2018: [0.0, 0.3555705652273914, 0.479317870754119, 0.5718070316694418, 1.0]
  Float_t binsx[5] = { 0.0, 0.36, 0.48, 0.57, 1.0 };
  histogram_mvaOutput_legacy_   = book1D(dir, "mvaOutput_legacy",   "mvaOutput_legacy",   4, binsx);

  histogram_mTauTauVis_ = book1D(dir, "mTauTauVis", "mTauTauVis", 20, 0., 200.);

  if(option_ == kOption_allHistograms)
  {
    histogram_memOutput_isValid_           = book1D(dir, "memOutput_isValid",           "memOutput_isValid",             3,  -1.5, +1.5);
    histogram_memOutput_errorFlag_         = book1D(dir, "memOutput_errorFlag",         "memOutput_errorFlag",           2,  -0.5, +1.5);
    histogram_memOutput_logWeight_ttH_     = book1D(dir, "memOutput_logWeight_ttH",     "memOutput_logWeight_ttH",     100, -20., +20.);
    histogram_memOutput_logWeight_ttZ_     = book1D(dir, "memOutput_logWeight_ttZ",     "memOutput_logWeight_ttZ",     100, -20., +20.);
    histogram_memOutput_logWeight_ttH_hww_ = book1D(dir, "memOutput_logWeight_ttH_hww", "memOutput_logWeight_ttH_hww", 100, -20., +20.);
    histogram_memOutput_LR_                = book1D(dir, "memOutput_LR",                "memOutput_LR",                 40,   0.,   1.);
    histogram_mem_logCPUTime_              = book1D(dir, "mem_logCPUTime",              "mem_logCPUTime",              400, -20., +20.);
    histogram_mem_logRealTime_             = book1D(dir, "mem_logRealTime",             "mem_logRealTime",             400, -20., +20.);
  }

  histogram_EventCounter_ = book1D(dir, "EventCounter", "EventCounter", 1, -0.5, +0.5);
}

void
EvtHistManager_3l_1tau::fillHistograms(const EvtHistManager_3l_1tau_Input & variables)
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
  }

  fillWithOverFlow(histogram_mvaOutput_legacy_, variables.mvaOutput_legacy, evtWeight, evtWeightErr);

  const double mTauTauVisSF = variables.mTauTauVis1 > 0. && variables.mTauTauVis2 > 0. ? 0.5 : 1.;
  if(variables.mTauTauVis1 > 0.)
  {
    fillWithOverFlow(histogram_mTauTauVis_, variables.mTauTauVis1, mTauTauVisSF * evtWeight, mTauTauVisSF * evtWeightErr);
  }
  if(variables.mTauTauVis2 > 0.)
  {
    fillWithOverFlow(histogram_mTauTauVis_, variables.mTauTauVis2, mTauTauVisSF * evtWeight, mTauTauVisSF * evtWeightErr);
  }

  if(option_ == kOption_allHistograms)
  {
    const MEMOutput_3l_1tau * const memOutput_3l_1tau = variables.memOutput_3l_1tau;
    if(memOutput_3l_1tau)
    {
      fillWithOverFlow(histogram_memOutput_isValid_, memOutput_3l_1tau->isValid(), evtWeight, evtWeightErr);

      if(memOutput_3l_1tau->isValid())
      {
        fillWithOverFlow(histogram_memOutput_errorFlag_, memOutput_3l_1tau->errorFlag(), evtWeight, evtWeightErr);

        if(memOutput_3l_1tau->errorFlag() == 0)
        {
          fillWithOverFlow(histogram_memOutput_logWeight_ttH_,     getLogWeight(memOutput_3l_1tau->weight_ttH()), evtWeight, evtWeightErr);
          fillWithOverFlow(histogram_memOutput_logWeight_ttZ_,     getLogWeight(memOutput_3l_1tau->weight_ttZ()), evtWeight, evtWeightErr);
          fillWithOverFlow(histogram_memOutput_logWeight_ttH_hww_, getLogWeight(memOutput_3l_1tau->weight_ttH_hww()), evtWeight, evtWeightErr);
          fillWithOverFlow(histogram_memOutput_LR_,                memOutput_3l_1tau->LR(), evtWeight, evtWeightErr);

          fillWithOverFlow(histogram_mem_logCPUTime_,  std::log(std::max(1.e-21f, memOutput_3l_1tau->cpuTime())),  evtWeight, evtWeightErr);
          fillWithOverFlow(histogram_mem_logRealTime_, std::log(std::max(1.e-21f, memOutput_3l_1tau->realTime())), evtWeight, evtWeightErr);
        }
      }
    }
    else
    {
      fillWithOverFlow(histogram_memOutput_isValid_, -1, evtWeight, evtWeightErr);
    }
  }

  fillWithOverFlow(histogram_EventCounter_, 0., evtWeight, evtWeightErr);
}
