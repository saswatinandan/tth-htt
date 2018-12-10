#ifndef tthAnalysis_HiggsToTauTau_HadTopTagger_semi_boosted_h
#define tthAnalysis_HiggsToTauTau_HadTopTagger_semi_boosted_h

#include "tthAnalysis/HiggsToTauTau/interface/RecoJetAK12.h"
#include "tthAnalysis/HiggsToTauTau/interface/RecoJet.h"
//#include "tthAnalysis/HiggsToTauTau/interface/hadTopTaggerAuxFunctions_geral.h" // kGenTop...
#include <map>

// forward declarations
//class HadTopKinFit;
class TMVAInterface;
class XGBInterface;

enum {
  kXGB_semi_boosted_with_kinFit, kXGB_semi_boosted_no_kinFit
};

class HadTopTagger_semi_boosted
{
public:
  HadTopTagger_semi_boosted(void) ;
  ~HadTopTagger_semi_boosted();

  /**
   * @brief Calculates MVA output.
   * @param mvaInputs Values of MVA input variables (stored in std::map with key = MVA input variable name)
   * @return          MVA output
   */
  std::map<int, double>
  operator()(
    const RecoJetAK12 & jet_ptrsAK12, const RecoJet & b_jet_candidate,
    bool & calculate_matching, bool & isGenMatched,
    double & genTopPt,
    std::map<int, Particle::LorentzVector> genVar, std::map<int, Particle::LorentzVector> genVarAnti
  );

  const std::map<std::string, double> &
  mvaInputs() const;

  //const HadTopKinFit *
  //kinFit() const;

protected:
  //HadTopKinFit * kinFit_;

  std::map<std::string, double> mvaInputsHTT;
  std::vector<std::string>      mvaInputsHTTSort;

  std::map<std::string, double> mvaInputsHTT_withKinFit;
  std::vector<std::string>      mvaInputsHTTSort_withKinFit;

  XGBInterface * mva_xgb_HTT_CSVsort4rd_;
  XGBInterface * mva_xgb_HTT_CSVsort4rd_withKinFit_;

};

#endif // tthAnalysis_HiggsToTauTau_HadTopTagger_h
