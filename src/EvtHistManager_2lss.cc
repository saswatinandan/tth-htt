#include "tthAnalysis/HiggsToTauTau/interface/EvtHistManager_2lss.h"

#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // fillWithOverFlow()
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // get_era(), kEra_*
#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h" // cmsException()

EvtHistManager_2lss::EvtHistManager_2lss(const edm::ParameterSet& cfg)
  : HistManagerBase(cfg)
  , era_(get_era(cfg.getParameter<std::string>("era")))
{}

const TH1 *
EvtHistManager_2lss::getHistogram_EventCounter() const
{
  return histogram_EventCounter_;
}

void EvtHistManager_2lss::bookHistograms(TFileDirectory & dir)
{
  histogram_numElectrons_    = book1D(dir, "numElectrons",    "numElectrons",     5, -0.5,  +4.5);
  histogram_numMuons_        = book1D(dir, "numMuons",        "numMuons",         5, -0.5,  +4.5);
  histogram_numHadTaus_      = book1D(dir, "numHadTaus",      "numHadTaus",       5, -0.5,  +4.5);
  histogram_numJets_         = book1D(dir, "numJets",         "numJets",         20, -0.5, +19.5);
  histogram_numBJets_loose_  = book1D(dir, "numBJets_loose",  "numBJets_loose",  10, -0.5,  +9.5);
  histogram_numBJets_medium_ = book1D(dir, "numBJets_medium", "numBJets_medium", 10, -0.5,  +9.5);
  histogram_numHTTv2_        = book1D(dir, "numHTTv2", "numHTTv2", 10, -0.5,  +9.5);

  histogram_numBJets_loose_vs_numJets_  = book2D(dir, "numBJets_loose_vs_numJets",  "numBJets_loose_vs_numJets",  8, -0.5, +7.5, 6, -0.5, +5.5);
  histogram_numBJets_medium_vs_numJets_ = book2D(dir, "numBJets_medium_vs_numJets", "numBJets_medium_vs_numJets", 8, -0.5, +7.5, 6, -0.5, +5.5);

  histogram_mvaOutput_2lss_ttV_   = book1D(dir, "mvaOutput_2lss_ttV",   "mvaOutput_2lss_ttV",   40, -1., +1.);
  histogram_output_NN_2lss_ttW_ttH_tH_3cat_v1_ = book1D(dir, "output_NN_2lss_ttW_ttH_tH_3cat_v1", "output_NN_2lss_ttW_ttH_tH_3cat_v1", 100, 0., +1.);
  histogram_mvaDiscr_2lss_        = book1D(dir, "mvaDiscr_2lss",        "mvaDiscr_2lss",         11,  -0.5, 10.5);

  histogram_output_NN_2lss_ttW_ttH_tH_3cat_v3_  = book1D(dir, "output_NN_2lss_ttW_ttH_tH_3cat_v3",  "output_NN_2lss_ttW_ttH_tH_3cat_v3",  100, 0, +1.);
  histogram_output_NN_2lss_ttW_ttH_3cat_v7_ = book1D(dir, "output_NN_2lss_ttW_ttH_3cat_v7", "output_NN_2lss_ttW_ttH_3cat_v7", 100, -1., +1.);

  histogram_EventCounter_ = book1D(dir, "EventCounter", "EventCounter", 1, -0.5, +0.5);

  histogram_output_NN_2lss_ttH_3cat_    = book1D(dir, "output_NN_2lss_ttH_3cat",    "output_NN_2lss_ttH_3cat",     100, 0.0,  +1.0);
  histogram_output_NN_2lss_ttW_ttH_3cat_ = book1D(dir, "output_NN_2lss_ttW_ttH_3cat",    "output_NN_2lss_ttW_ttH_3cat",     100, 0.0,  +1.0);

  histogram_output_NN_2lss_ttW_ttH_tH_4cat_v1_    = book1D(dir, "output_NN_2lss_ttW_ttH_tH_4cat_v1",    "output_NN_2lss_ttW_ttH_tH_4cat_v1",     100, 0.0,  +1.0);
  histogram_output_NN_2lss_ttW_ttH_tH_4cat_v2_    = book1D(dir, "output_NN_2lss_ttW_ttH_tH_4cat_v2",    "output_NN_2lss_ttW_ttH_tH_4cat_v2",     100, 0.0,  +1.0);
  histogram_mva_Updated_    = book1D(dir, "mva_Updated",    "mva_Updated",     100, 0.0,  +1.0);
  histogram_mva_oldVar_    = book1D(dir, "mva_oldVar",    "mva_oldVar",     100, 0.0,  +1.0);

}

void
EvtHistManager_2lss::fillHistograms(int numElectrons,
                                    int numMuons,
                                    int numHadTaus,
                                    int numJets,
                                    int numBJets_loose,
                                    int numBJets_medium,
                                    int numHTTv2,
                                    double evtWeight,
                                    double mvaOutput_2lss_ttV,
                                    double output_NN_2lss_ttW_ttH_tH_3cat_v1,
                                    double mvaDiscr_2lss,
                                    double output_NN_2lss_ttW_ttH_tH_3cat_v3,
                                    double output_NN_2lss_ttW_ttH_3cat_v7,
                                    //
                                    double output_NN_2lss_ttH_3cat,
                                    double output_NN_2lss_ttW_ttH_3cat,
                                    double output_NN_2lss_ttW_ttH_tH_4cat_v1,
                                    double output_NN_2lss_ttW_ttH_tH_4cat_v2,
                                    double mva_Updated,
                                    double mva_oldVar
                                  )
{
  const double evtWeightErr = 0.;

  fillWithOverFlow(histogram_numElectrons_,    numElectrons,    evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_numMuons_,        numMuons,        evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_numHadTaus_,      numHadTaus,      evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_numJets_,         numJets,         evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_numBJets_loose_,  numBJets_loose,  evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_numBJets_medium_, numBJets_medium, evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_numHTTv2_,        numHTTv2, evtWeight, evtWeightErr);

  fillWithOverFlow2d(histogram_numBJets_loose_vs_numJets_,  numJets, numBJets_loose,  evtWeight, evtWeightErr);
  fillWithOverFlow2d(histogram_numBJets_medium_vs_numJets_, numJets, numBJets_medium, evtWeight, evtWeightErr);

  fillWithOverFlow(histogram_mvaOutput_2lss_ttV_,   mvaOutput_2lss_ttV,   evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_output_NN_2lss_ttW_ttH_tH_3cat_v1_, output_NN_2lss_ttW_ttH_tH_3cat_v1, evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_mvaDiscr_2lss_,        mvaDiscr_2lss,        evtWeight, evtWeightErr);

  fillWithOverFlow(histogram_output_NN_2lss_ttW_ttH_tH_3cat_v3_,  output_NN_2lss_ttW_ttH_tH_3cat_v3,  evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_output_NN_2lss_ttW_ttH_3cat_v7_, output_NN_2lss_ttW_ttH_3cat_v7, evtWeight, evtWeightErr);

  fillWithOverFlow(histogram_EventCounter_, 0., evtWeight, evtWeightErr);

  fillWithOverFlow(histogram_output_NN_2lss_ttH_3cat_,           output_NN_2lss_ttH_3cat,    evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_output_NN_2lss_ttW_ttH_3cat_,   output_NN_2lss_ttW_ttH_3cat,        evtWeight, evtWeightErr);

  fillWithOverFlow(histogram_output_NN_2lss_ttW_ttH_tH_4cat_v1_,  output_NN_2lss_ttW_ttH_tH_4cat_v1,  evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_output_NN_2lss_ttW_ttH_tH_4cat_v2_, output_NN_2lss_ttW_ttH_tH_4cat_v2, evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_mva_Updated_,        mva_Updated, evtWeight, evtWeightErr);
  fillWithOverFlow(histogram_mva_oldVar_,        mva_oldVar, evtWeight, evtWeightErr);

}
