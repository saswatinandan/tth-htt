#ifndef tthAnalysis_HiggsToTauTau_MEMOutputReader_3l_1tau_h
#define tthAnalysis_HiggsToTauTau_MEMOutputReader_3l_1tau_h

#include "tthAnalysis/HiggsToTauTau/interface/ReaderBase.h" // ReaderBase
#include "tthAnalysis/HiggsToTauTau/interface/MEMOutput_3l_1tau.h" // MEMOutput_3l_1tau

// forward declarations
class TTree;

class MEMOutputReader_3l_1tau
  : public ReaderBase
{
public:
  MEMOutputReader_3l_1tau(const std::string & branchName_num,
                          const std::string & branchName_obj);
  ~MEMOutputReader_3l_1tau();

  /**
   * @brief Call tree->SetBranchAddress for all MEMOutput_3l_1tau branches
   */
  std::vector<std::string>
  setBranchAddresses(TTree * tree) override;

  /**
   * @brief Read branches from tree and use information to fill collection of MEMOutput_3l_1tau objects
   * @return Collection of MEMOutput_3l_1tau objects
   */
  std::vector<MEMOutput_3l_1tau>
  read() const;

protected:
 /**
   * @brief Initialize names of branches to be read from tree
   */
  void
  setBranchNames();

  const int max_nMEMOutputs_;
  std::string branchName_num_;
  std::string branchName_obj_;

  std::string branchName_run_;
  std::string branchName_lumi_;
  std::string branchName_evt_;
  std::string branchName_leadLepton_eta_;
  std::string branchName_leadLepton_phi_;
  std::string branchName_subleadLepton_eta_;
  std::string branchName_subleadLepton_phi_;
  std::string branchName_thirdLepton_eta_;
  std::string branchName_thirdLepton_phi_;
  std::string branchName_hadTau_eta_;
  std::string branchName_hadTau_phi_;
  std::string branchName_weight_ttH_;
  std::string branchName_weight_ttH_error_;
  std::string branchName_weight_ttZ_;
  std::string branchName_weight_ttZ_error_;
  std::string branchName_weight_ttH_hww_;
  std::string branchName_weight_ttH_hww_error_;
  std::string branchName_LR_;
  std::string branchName_LR_up_;
  std::string branchName_LR_down_;
  std::string branchName_cpuTime_;
  std::string branchName_realTime_;
  std::string branchName_isValid_;
  std::string branchName_errorFlag_;

  Int_t nMEMOutputs_;
  UInt_t * run_;
  UInt_t * lumi_;
  ULong64_t * evt_;
  Float_t * leadLepton_eta_;
  Float_t * leadLepton_phi_;
  Float_t * subleadLepton_eta_;
  Float_t * subleadLepton_phi_;
  Float_t * thirdLepton_eta_;
  Float_t * thirdLepton_phi_;
  Float_t * hadTau_eta_;
  Float_t * hadTau_phi_;
  Float_t * weight_ttH_;
  Float_t * weight_ttH_error_;
  Float_t * weight_ttZ_;
  Float_t * weight_ttZ_error_;
  Float_t * weight_ttH_hww_;
  Float_t * weight_ttH_hww_error_;
  Float_t * LR_;
  Float_t * LR_up_;
  Float_t * LR_down_;
  Float_t * cpuTime_;
  Float_t * realTime_;
  Int_t * isValid_;
  Int_t * errorFlag_;

  // CV: make sure that only one MEMOutputReader_3l_1tau instance exists for a given branchName,
  //     as ROOT cannot handle multiple TTree::SetBranchAddress calls for the same branch.
  static std::map<std::string, int> numInstances_;
  static std::map<std::string, MEMOutputReader_3l_1tau *> instances_;
};

#endif // tthAnalysis_HiggsToTauTau_MEMOutputReader_3l_1tau_h

