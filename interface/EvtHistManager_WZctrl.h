#ifndef tthAnalysis_HiggsToTauTau_EvtHistManager_WZctrl_h
#define tthAnalysis_HiggsToTauTau_EvtHistManager_WZctrl_h

/** \class EvtHistManager_WZctrl
 *
 * Book and fill histograms for event-level quantities in ttH, H->tautau analysis
 * in control region enriched in WZ background
 *
 * \author Christian Veelken, Tallinn
 *
 */

#include "tthAnalysis/HiggsToTauTau/interface/HistManagerBase.h"      // HistManagerBase
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // Era

class EvtHistManager_WZctrl
  : public HistManagerBase
{
public:
  EvtHistManager_WZctrl(const edm::ParameterSet & cfg);
  ~EvtHistManager_WZctrl() {}

  /// book and fill histograms
  void
  bookHistograms(TFileDirectory & dir) override;

  void
  fillHistograms(int numElectrons,
                 int numMuons,
                 int numHadTaus,
                 int numJets,
                 int numBJets_loose,
                 int numBJets_medium,
                 double mvaOutput_2lss_ttV,
                 double mvaOutput_2lss_ttbar,
                 double mvaDiscr_2lss,
                 double mvaOutput_3l_ttV,
                 double mvaOutput_3l_ttbar,
                 double mvaDiscr_3l,
                 double mLL,
                 double mT,
                 int sumLeptonCharge,
		 double met,
                 double evtWeight);

private:
  TH1 * histogram_numElectrons_;
  TH1 * histogram_numMuons_;
  TH1 * histogram_numLeptons_;
  TH1 * histogram_numHadTaus_;
  TH1 * histogram_numJets_;
  TH1 * histogram_numBJets_loose_;
  TH1 * histogram_numBJets_medium_;

  TH1 * histogram_mvaOutput_2lss_ttV_;
  TH1 * histogram_mvaOutput_2lss_ttbar_;
  TH1 * histogram_mvaDiscr_2lss_;

  TH1 * histogram_mvaOutput_3l_ttV_;
  TH1 * histogram_mvaOutput_3l_ttbar_;
  TH1 * histogram_mvaDiscr_3l_;

  TH1 * histogram_mLL_;
  TH1 * histogram_mT_;
  TH1 * histogram_sumLeptonCharge_;
  TH1 * histogram_met_;
  
  TH1 * histogram_EventCounter_;
};

#endif
