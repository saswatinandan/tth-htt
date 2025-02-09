#include "tthAnalysis/HiggsToTauTau/interface/EvtHistManager_3l.h"

#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // fillWithOverFlow(), getLogWeight()
#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h" // cmsException()

EvtHistManager_3l::EvtHistManager_3l(const edm::ParameterSet & cfg, bool isControlRegion)
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
    "mvaOutput_3l_ttV",
    "mvaOutput_3l_ttbar",
    "memOutput_isValid",
    "memOutput_errorFlag",
    "memOutput_logWeight_ttH",
    "memOutput_logWeight_tt",
    "memOutput_LR",
    "mem_logCPUTime",
    "mem_logRealTime",
    "output_NN_cr"
  };
  std::vector<std::string> sysOpts_all = {
    "EventCounter"
  };
  if(isControlRegion)
  {
    sysOpts_all.push_back("control_eee");
    sysOpts_all.push_back("control_eem");
    sysOpts_all.push_back("control_emm");
    sysOpts_all.push_back("control_mmm");
  }
  else
  {
    sysOpts_all.push_back("output_NN_ttH_bl");
    sysOpts_all.push_back("output_NN_ttH_bt");
    sysOpts_all.push_back("output_NN_tH_bl");
    sysOpts_all.push_back("output_NN_tH_bt");
    sysOpts_all.push_back("output_NN_rest_eee");
    sysOpts_all.push_back("output_NN_rest_eem");
    sysOpts_all.push_back("output_NN_rest_emm");
    sysOpts_all.push_back("output_NN_rest_mmm");

    sysOpts_all.push_back("mass_3L_lj_pos");
    sysOpts_all.push_back("mass_3L_lj_neg");
    sysOpts_all.push_back("mass_3L_hj_pos");
    sysOpts_all.push_back("mass_3L_hj_neg");

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
EvtHistManager_3l::getHistogram_EventCounter() const
{
  return histogram_EventCounter_;
}

void
EvtHistManager_3l::bookCategories(TFileDirectory & dir,
                                  const std::map<std::string, std::vector<double>> & categories_list_NN,
                                  const std::map<std::string, std::vector<double>> & categories_list_SVA,
                                  bool isControlRegion)
{
  for(auto category: categories_list_NN)
  {
    if(isControlRegion)
    {
      histograms_by_category_[category.first] = book1D(dir, category.first, category.first, 12,  0.5, +12.5);
    }
    else if(! category.second.empty())
    {
      const int npoints = category.second.size();
      Float_t binsx[npoints];
      std::copy(category.second.begin(), category.second.end(), binsx);
      histograms_by_category_[category.first] = book1D(dir, category.first, category.first, npoints - 1, binsx);
    }
    else
    {
      histograms_by_category_[category.first] = book1D(dir, category.first, category.first, 100,  0., +1.);
    }
    central_or_shiftOptions_[category.first] = { "*" };
  }
  if(! isControlRegion)
  {
    Float_t bins_mass_3L[6] = {20.,100.,140.,190.,250.,1000.};
    for(auto category: categories_list_SVA)
    {
      histograms_by_category_SVA_[category.first] = book1D(dir, category.first, category.first, 5,  bins_mass_3L);
      central_or_shiftOptions_[category.first] = { "*" };
    }
  }
}

void
EvtHistManager_3l::bookHistograms(TFileDirectory & dir)
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

  histogram_mvaOutput_3l_ttV_   = book1D(dir, "mvaOutput_3l_ttV",   "mvaOutput_3l_ttV",   40, -1., +1.);
  histogram_mvaOutput_3l_ttbar_ = book1D(dir, "mvaOutput_3l_ttbar", "mvaOutput_3l_ttbar", 40, -1., +1.);

  if(option_ == kOption_allHistograms)
  {
    histogram_memOutput_isValid_       = book1D(dir, "memOutput_isValid",       "memOutput_isValid",         3,  -1.5, +1.5);
    histogram_memOutput_errorFlag_     = book1D(dir, "memOutput_errorFlag",     "memOutput_errorFlag",       2,  -0.5, +1.5);
    histogram_memOutput_logWeight_ttH_ = book1D(dir, "memOutput_logWeight_ttH", "memOutput_logWeight_ttH", 100, -20., +20.);
    histogram_memOutput_logWeight_tt_  = book1D(dir, "memOutput_logWeight_tt",  "memOutput_logWeight_tt",  100, -20., +20.);
    histogram_memOutput_LR_            = book1D(dir, "memOutput_LR",            "memOutput_LR",             40,   0.,   1.);
    histogram_mem_logCPUTime_          = book1D(dir, "mem_logCPUTime",          "mem_logCPUTime",          400, -20., +20.);
    histogram_mem_logRealTime_         = book1D(dir, "mem_logRealTime",         "mem_logRealTime",         400, -20., +20.);
  }

  histogram_EventCounter_ = book1D(dir, "EventCounter", "EventCounter", 1, -0.5, +0.5);
}

void
EvtHistManager_3l::fillHistograms(const EvtHistManager_3l_Input & variables)
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

  fillWithOverFlow(histogram_mvaOutput_3l_ttV_,   variables.mvaOutput_3l_ttV,   evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_mvaOutput_3l_ttbar_, variables.mvaOutput_3l_ttbar, evtWeight, evtWeightErr);

  if(! histograms_by_category_.count(variables.category_NN))
  {
    throw cmsException(this, __func__, __LINE__) << "Histogram of the name '" << variables.category_NN << "' was never booked";
  }
  fillWithOverFlow(histograms_by_category_[variables.category_NN], variables.mvaOutput_category_NN, evtWeight, evtWeightErr);

  if(! variables.isControlRegion)
  {
    if(! histograms_by_category_SVA_.count(variables.category_SVA))
    {
      throw cmsException(this, __func__, __LINE__) << "Histogram of the name '" << variables.category_SVA << "' was never booked";
    }
    fillWithOverFlow(histograms_by_category_SVA_[variables.category_SVA], variables.mass_3L, evtWeight, evtWeightErr);
  }

  if(option_ == kOption_allHistograms)
  {
    const MEMOutput_3l * const memOutput_3l = variables.memOutput_3l;
    if(memOutput_3l)
    {
      fillWithOverFlow(histogram_memOutput_isValid_, memOutput_3l->isValid(), evtWeight, evtWeightErr);

      if(memOutput_3l->isValid())
      {
        fillWithOverFlow(histogram_memOutput_errorFlag_, memOutput_3l->errorFlag(), evtWeight, evtWeightErr);

        if(memOutput_3l->errorFlag() == 0)
        {
          fillWithOverFlow(histogram_memOutput_logWeight_ttH_, getLogWeight(memOutput_3l->weight_ttH()), evtWeight, evtWeightErr);
          fillWithOverFlow(histogram_memOutput_logWeight_tt_,  getLogWeight(memOutput_3l->weight_tt()),  evtWeight, evtWeightErr);
          fillWithOverFlow(histogram_memOutput_LR_,            memOutput_3l->LR(),                       evtWeight, evtWeightErr);

          fillWithOverFlow(histogram_mem_logCPUTime_,  std::log(std::max(1.e-21f, memOutput_3l->cpuTime())),  evtWeight, evtWeightErr);
          fillWithOverFlow(histogram_mem_logRealTime_, std::log(std::max(1.e-21f, memOutput_3l->realTime())), evtWeight, evtWeightErr);
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
