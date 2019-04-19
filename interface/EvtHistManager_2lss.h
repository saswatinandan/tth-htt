#ifndef tthAnalysis_HiggsToTauTau_EvtHistManager_2lss_h
#define tthAnalysis_HiggsToTauTau_EvtHistManager_2lss_h

/** \class EvtHistManager_2lss
 *
 * Book and fill histograms for event-level quantities in ttH multilepton analysis
 * in 2lss category
 *
 * \author Christian Veelken, Tallin
 *
 */

#include "tthAnalysis/HiggsToTauTau/interface/HistManagerBase.h" // HistManagerBase

class EvtHistManager_2lss
  : public HistManagerBase
{
 public:
  EvtHistManager_2lss(const edm::ParameterSet & cfg);
  ~EvtHistManager_2lss() {}

  void
  bookHistograms(TFileDirectory & dir) override;

  void
  fillHistograms(int numElectrons,
                 int numMuons,
                 int numHadTaus,
                 int numJets,
                 int numBJets_loose,
                 int numBJets_medium,
                 int num_HTTv2,
                 double evtWeight,
                 //
                 double mvaOutput_2lss_ttV,
                 double output_NN_2lss_ttH_tH_4cat_onlyTHQ_v5,
                 double mvaDiscr_2lss,
                 double output_NN_2lss_ttH_tH_4cat_onlyTHQ_v7,
                 double output_NN_2lss_ttH_tH_4cat_onlyTHQ_v3,
                 //
                 double output_NN_2lss_ttH_tH_4cat_onlyTHQ_v1,
                 double output_NN_2lss_ttH_tH_4cat_onlyTHQ_v2,
                 double output_NN_2lss_ttH_tH_4cat_onlyTHQ_v4,
                 double output_NN_2lss_ttH_tH_4cat_onlyTHQ_v6,
                 double mva_Updated,
                 double mva_oldVar
               );

  const TH1 *
  getHistogram_EventCounter() const;

 private:
  int era_;

  TH1 * histogram_numElectrons_;
  TH1 * histogram_numMuons_;
  TH1 * histogram_numHadTaus_;
  TH1 * histogram_numJets_;
  TH1 * histogram_numBJets_loose_;
  TH1 * histogram_numBJets_medium_;
  TH1 * histogram_numHTTv2_;

  // CV: used to check loss in signal efficiency in case events with
  // high jet and b-jet multiplicity are vetoed to avoid overlap with ttH, H->bb analysis
  // (alternative: ttH, H->bb analysis adds hadronic tau veto)
  TH2 * histogram_numBJets_loose_vs_numJets_;
  TH2 * histogram_numBJets_medium_vs_numJets_;

  TH1 * histogram_mvaOutput_2lss_ttV_;
  TH1 * histogram_mvaDiscr_2lss_;

  TH1 * histogram_output_NN_2lss_ttH_tH_4cat_onlyTHQ_v7_;

  TH1 * histogram_EventCounter_;

  TH1 * histogram_output_NN_2lss_ttH_tH_4cat_onlyTHQ_v1_;
  TH1 * histogram_output_NN_2lss_ttH_tH_4cat_onlyTHQ_v2_;
  TH1 * histogram_output_NN_2lss_ttH_tH_4cat_onlyTHQ_v3_;
  TH1 * histogram_output_NN_2lss_ttH_tH_4cat_onlyTHQ_v4_;
  TH1 * histogram_output_NN_2lss_ttH_tH_4cat_onlyTHQ_v5_;
  TH1 * histogram_output_NN_2lss_ttH_tH_4cat_onlyTHQ_v6_;
  TH1 * histogram_mva_Updated_;
  TH1 * histogram_mva_oldVar_;

};

#endif
