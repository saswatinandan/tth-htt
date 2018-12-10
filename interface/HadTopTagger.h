#ifndef tthAnalysis_HiggsToTauTau_HadTopTagger_h
#define tthAnalysis_HiggsToTauTau_HadTopTagger_h

#include "tthAnalysis/HiggsToTauTau/interface/RecoJet.h" // RecoJet

// forward declarations
class HadTopKinFit;
class TMVAInterface;
class XGBInterface;

enum {
  kXGB_with_kinFit, kXGB_multilep,
  kXGB_CSVsort3rd, kXGB_CSVsort3rd_withKinFit,
  kXGB_highestCSV, kXGB_highestCSV_withKinFit
};

class HadTopTagger
{
public:
  HadTopTagger(void) ;
  ~HadTopTagger();

  /**
   * @brief Calculates MVA output.
   * @param mvaInputs Values of MVA input variables (stored in std::map with key = MVA input variable name)
   * @return          MVA output
   */
  std::map<int, double>
  operator()(const RecoJet & recBJet,
             const RecoJet & recWJet1,
             const RecoJet & recWJet2,
             bool & calculate_matching, bool & isGenMatched,
             double & genTopPt,
             std::map<int, Particle::LorentzVector> genVar, std::map<int, Particle::LorentzVector> genVarAnti,
             bool massCut
           );

  const std::map<std::string, double> &
  mvaInputs() const;

protected:
  HadTopKinFit * kinFit_;

  std::map<std::string, double> mvaInputsWithKinFit;
  std::vector<std::string>      mvaInputsWithKinFitSort;
  TMVAInterface * mva_hadTopTagger_xgb_withKinFit_;

  std::map<std::string, double> mvaInputsHTT_multilep;
  std::vector<std::string>      mvaInputsHTT_multilepSort;
  TMVAInterface * mva_hadTopTagger_xgb_HTT_multilep_;

  std::map<std::string, double> mvaInputsHTT;
  std::vector<std::string>      mvaInputsHTTSort;

  std::map<std::string, double> mvaInputsHTT_withKinFit;
  std::vector<std::string>      mvaInputsHTT_withKinFitSort;

  XGBInterface * mva_xgb_HTT_CSVsort3rd_;
  XGBInterface * mva_xgb_HTT_CSVsort3rd_withKinFit_;
  XGBInterface * mva_xgb_HTT_highestCSV_;
  XGBInterface * mva_xgb_HTT_highestCSV_withKinFit_;

};

#endif // tthAnalysis_HiggsToTauTau_HadTopTagger_h
