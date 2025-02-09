#ifndef tthAnalysis_HiggsToTauTau_MEtFilterReader_h
#define tthAnalysis_HiggsToTauTau_MEtFilterReader_h

#include "tthAnalysis/HiggsToTauTau/interface/ReaderBase.h"           // ReaderBase
#include "tthAnalysis/HiggsToTauTau/interface/MEtFilterFlag.h"        // MEtFilterFlag::
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // Era

#include <array> // std::array<>

// forward declarations
class MEtFilter;

class MEtFilterReader
  : public ReaderBase
{
public:
  MEtFilterReader(MEtFilter * metFilter,
                  Era era);
  ~MEtFilterReader();

  std::vector<std::string>
  setBranchAddresses(TTree * tree) override;

protected:
 /**
   * @brief Initialize names of branches to be read from tree
   */
  void
  setBranchNames();

  std::array<std::string, MEtFilterFlag::LAST> branchNames_;

  // CV: make sure that only one MEtFilterReader instance exists,
  static int numInstances_;
  static MEtFilterReader * instance_;

  MEtFilter * metFilter_;
  Era era_;
};

#endif
