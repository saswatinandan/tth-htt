#ifndef tthAnalysis_HiggsToTauTau_EvtHistManager_3l_h
#define tthAnalysis_HiggsToTauTau_EvtHistManager_3l_h

/** \class EvtHistManager_3l
 *
 * Book and fill histograms for event-level quantities in ttH multilepton analysis
 * in 3l category
 *
 * \author Christian Veelken, Tallin
 *
 */

#include "tthAnalysis/HiggsToTauTau/interface/HistManagerBase.h" // HistManagerBase

#include "tthAnalysis/HiggsToTauTau/interface/MEMOutput_3l.h" // MEMOutput_3l

struct EvtHistManager_3l_Input
{
  std::size_t numElectrons;
  std::size_t numMuons;
  std::size_t numHadTaus;
  std::size_t numJets;
  std::size_t numBJets_loose;
  std::size_t numBJets_medium;
  double mvaOutput_3l_ttV;
  double mvaOutput_3l_ttbar;
  double mass_3L;
  std::string category_SVA;
  double mvaOutput_category_NN;
  std::string category_NN;
  const MEMOutput_3l * memOutput_3l;
  double evtWeight;
  bool isControlRegion;
};

class EvtHistManager_3l
  : public HistManagerBase
{
 public:
  EvtHistManager_3l(const edm::ParameterSet & cfg, bool isControlRegion);
  ~EvtHistManager_3l() {}

  /// book and fill histograms
  void
  bookHistograms(TFileDirectory & dir) override;

  void
  bookCategories(TFileDirectory & dir,
                 const std::map<std::string, std::vector<double>> & categories_list_NN,
                 const std::map<std::string, std::vector<double>> & categories_list_SVA,
                 bool isControlRegion);

  void
  fillHistograms(const EvtHistManager_3l_Input & variables);

  const TH1 *
  getHistogram_EventCounter() const;

  enum { kOption_undefined, kOption_allHistograms, kOption_minimalHistograms };

 private:
  std::vector<std::string> ctrl_cateories_;

  TH1 * histogram_numElectrons_;
  TH1 * histogram_numMuons_;
  TH1 * histogram_numHadTaus_;
  TH1 * histogram_numJets_;
  TH1 * histogram_numBJets_loose_;
  TH1 * histogram_numBJets_medium_;

  // CV: used to check loss in signal efficiency in case events with high jet and b-jet multiplicity are vetoed
  // to avoid overlap with ttH, H->bb analysis (alternative: ttH, H->bb analysis adds hadronic tau veto)
  TH2 * histogram_numBJets_loose_vs_numJets_;
  TH2 * histogram_numBJets_medium_vs_numJets_;

  TH1 * histogram_mvaOutput_3l_ttV_;
  TH1 * histogram_mvaOutput_3l_ttbar_;

  std::map<std::string, TH1 *> histograms_by_category_;
  std::map<std::string, TH1 *> histograms_by_category_SVA_;

  TH1 * histogram_memOutput_isValid_;
  TH1 * histogram_memOutput_errorFlag_;
  TH1 * histogram_memOutput_logWeight_ttH_;
  TH1 * histogram_memOutput_logWeight_tt_;
  TH1 * histogram_memOutput_LR_;
  TH1 * histogram_mem_logCPUTime_;
  TH1 * histogram_mem_logRealTime_;

  TH1 * histogram_EventCounter_;
  int option_; // flag to book & fill either full or minimal set of histograms (to reduce memory consumption of hadd jobs)
};

#endif
